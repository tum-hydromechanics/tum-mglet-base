#!/usr/bin/env python3
"""Local correctness and performance harness for MGLET particles."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shutil
import socket
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

try:
    import h5py
    import numpy as np
except ImportError:
    h5py = None
    np = None


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
WORKSPACE_ROOT = REPO_ROOT.parent
DEFAULT_MANIFEST = SCRIPT_DIR / "benchmarks.json"
DEFAULT_RESULTS_ROOT = WORKSPACE_ROOT / "benchmark-results"
BCC_EXTERNAL_ROOT = Path("/home/yaydin/thesis/particle_run/10_BCC")
OUTPUT_NAMES = {
    "mglet-perf-report.txt",
    "mglet.OUT",
    "particles.h5",
    "fields.h5",
    "ib_stencils.h5",
}
OUTPUT_DIRS = {"LOGS", "Particle_Snapshots", "Particle_Statistics"}


class BenchmarkError(RuntimeError):
    """Raised when a correctness or benchmark run is invalid."""


def require_within_workspace(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(WORKSPACE_ROOT)
    except ValueError as exc:
        raise BenchmarkError(
            f"{label} must be inside {WORKSPACE_ROOT}, got {resolved}"
        ) from exc
    return resolved


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise BenchmarkError(f"Expected a JSON object in {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def deep_merge(target: dict[str, Any], updates: dict[str, Any]) -> None:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            deep_merge(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def load_manifest(path: Path) -> dict[str, Any]:
    manifest = read_json(path)
    if manifest.get("version") != 1:
        raise BenchmarkError(f"Unsupported manifest version in {path}")
    for key in ("timers", "correctness", "benchmarks"):
        if key not in manifest or not isinstance(manifest[key], list):
            raise BenchmarkError(f"Manifest key {key!r} must be a list")
        if key != "timers":
            for entry in manifest[key]:
                if not isinstance(entry, dict):
                    raise BenchmarkError(f"Manifest {key} entries must be objects")
                safe_name(entry.get("name"), f"{key} case name")
    return manifest


def safe_name(value: Any, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value in {".", ".."}
        or Path(value).name != value
        or "/" in value
        or "\\" in value
    ):
        raise BenchmarkError(f"{label} must be one safe path component: {value!r}")
    return value


def source_path(entry: dict[str, Any]) -> Path:
    if "external_root" in entry:
        root = Path(entry["external_root"]).expanduser().resolve()
        if root != BCC_EXTERNAL_ROOT:
            raise BenchmarkError(f"External root is not permitted: {root}")
        source = (root / entry["source"]).resolve()
        permitted_root = root
    else:
        source = (SCRIPT_DIR / entry["source"]).resolve()
        permitted_root = REPO_ROOT
    try:
        source.relative_to(permitted_root)
    except ValueError as exc:
        raise BenchmarkError(f"Case source escapes permitted root: {source}") from exc
    if not (source / "parameters.json").is_file():
        raise BenchmarkError(f"Case has no parameters.json: {source}")
    return source


def source_identifier(entry: dict[str, Any]) -> str:
    source = source_path(entry)
    if "external_root" in entry:
        return f"10_BCC/{source.relative_to(BCC_EXTERNAL_ROOT)}"
    return str(source.relative_to(REPO_ROOT))


def external_input_paths(entry: dict[str, Any]) -> dict[str, Path]:
    bindings = entry.get("external_inputs", {})
    if not isinstance(bindings, dict):
        raise BenchmarkError("external_inputs must be an object")
    resolved = {}
    for alias, relative_path in bindings.items():
        safe_name(alias, "external input alias")
        if not isinstance(relative_path, str):
            raise BenchmarkError("External input paths must be strings")
        path = (BCC_EXTERNAL_ROOT / relative_path).resolve()
        try:
            path.relative_to(BCC_EXTERNAL_ROOT)
        except ValueError as exc:
            raise BenchmarkError(f"External input escapes 10_BCC: {path}") from exc
        if not path.is_file():
            raise BenchmarkError(f"External input does not exist: {path}")
        resolved[alias] = path
    return resolved


def selected_h5_inputs(parameters: dict[str, Any]) -> set[str]:
    selected: set[str] = set()
    io_config = parameters.get("io", {})
    grid_file = io_config.get("grids")
    if isinstance(grid_file, str):
        selected.add(grid_file)
    time_config = parameters.get("time", {})
    if time_config.get("read", False):
        infile = io_config.get("infile")
        if isinstance(infile, str):
            selected.add(infile)
    ib_config = parameters.get("ib", {})
    if ib_config.get("type", "noib") != "noib":
        stencil_file = ib_config.get("stencilfile", "ib_stencils.h5")
        if isinstance(stencil_file, str):
            selected.add(stencil_file)
    return selected


def validate_parameter_paths(parameters: dict[str, Any]) -> None:
    io_config = parameters.get("io", {})
    if not isinstance(io_config, dict):
        raise BenchmarkError("Parameter section 'io' must be an object")
    for key in ("grids", "infile", "outfile"):
        value = io_config.get(key)
        if value is not None:
            safe_name(value, f"io.{key}")
    ib_config = parameters.get("ib", {})
    if not isinstance(ib_config, dict):
        raise BenchmarkError("Parameter section 'ib' must be an object")
    stencil_file = ib_config.get("stencilfile")
    if stencil_file is not None:
        safe_name(stencil_file, "ib.stencilfile")
    for geometry in ib_config.get("geometries", []):
        if not isinstance(geometry, dict):
            raise BenchmarkError("Each IB geometry must be an object")
        geometry_file = geometry.get("file")
        if geometry_file is not None:
            safe_name(geometry_file, "ib.geometries.file")


def apply_generated_seed(entry: dict[str, Any], parameters: dict[str, Any]) -> None:
    if "seed_base" not in entry:
        return
    words_per_rank = int(entry.get("seed_words_per_rank", 8))
    ranks = int(entry["ranks"])
    base = int(entry["seed_base"])
    parameters.setdefault("particles", {})["particle_seed"] = [
        base + 104729 * index for index in range(words_per_rank * ranks)
    ]


def copy_case_input(source: Path, destination: Path, entry: dict[str, Any]) -> None:
    if source.name != "ParticleDict.txt" or entry.get(
        "particle_dict_columns", "id_xyz"
    ) == "id_xyz":
        shutil.copy2(source, destination)
        return
    if entry["particle_dict_columns"] != "xyz":
        raise BenchmarkError("particle_dict_columns must be 'id_xyz' or 'xyz'")
    with source.open(encoding="utf-8") as input_handle, destination.open(
        "w", encoding="utf-8"
    ) as output_handle:
        count_line = input_handle.readline()
        if not count_line:
            raise BenchmarkError(f"Empty particle dictionary: {source}")
        output_handle.write(count_line)
        for particle_id, line in enumerate(input_handle, start=1):
            output_handle.write(f"{particle_id} {line}")


def stage_case(entry: dict[str, Any], destination: Path) -> dict[str, Any]:
    source = source_path(entry)
    parameters = read_json(source / "parameters.json")
    for key in entry.get("remove_top_level", []):
        parameters.pop(key, None)
    deep_merge(parameters, entry.get("overrides", {}))
    apply_generated_seed(entry, parameters)
    validate_parameter_paths(parameters)

    destination.mkdir(parents=True, exist_ok=False)
    selected_h5 = selected_h5_inputs(parameters)
    copy_inputs = entry.get("copy_inputs")
    if copy_inputs is None:
        source_items = list(source.iterdir())
    else:
        source_items = []
        for name in copy_inputs:
            safe_name(name, "copied input")
            item = source / name
            if not item.is_file():
                raise BenchmarkError(f"Copied input does not exist: {item}")
            source_items.append(item)
    for item in source_items:
        if item.is_file() and item.name != "parameters.json":
            if item.suffix != ".h5" or item.name in selected_h5:
                copy_case_input(item, destination / item.name, entry)
    for alias, target in external_input_paths(entry).items():
        (destination / alias).symlink_to(target)
    for name in selected_h5:
        if not (destination / name).is_file():
            raise BenchmarkError(f"Staged HDF5 input is missing: {name}")
    write_json(destination / "parameters.json", parameters)
    return parameters


def command_output(command: list[str], cwd: Path | None = None) -> str:
    try:
        return subprocess.run(
            command,
            cwd=cwd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def run_mglet(
    binary: Path,
    run_dir: Path,
    ranks: int,
    mpirun: str,
    mpirun_args: list[str],
) -> dict[str, Any]:
    command = [mpirun, *mpirun_args, "-n", str(ranks), str(binary), "parameters.json"]
    log_path = run_dir / "mglet.OUT"
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"COMMAND: {' '.join(command)}\n\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=run_dir,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    wall_seconds = time.perf_counter() - started
    if completed.returncode != 0:
        raise BenchmarkError(
            f"MGLET failed with status {completed.returncode}; see {log_path}"
        )
    return {
        "command": command,
        "returncode": completed.returncode,
        "wall_seconds": wall_seconds,
        "output_bytes": generated_output_size(run_dir),
    }


def generated_output_size(run_dir: Path) -> int:
    total = 0
    for path in run_dir.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(run_dir)
        if relative.parts[0] in OUTPUT_DIRS or path.name in OUTPUT_NAMES:
            total += path.stat().st_size
    return total


def directory_size(path: Path) -> int:
    if not path.is_dir():
        return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def parse_timer_report(path: Path, timer_ids: set[int]) -> dict[str, Any]:
    if not path.is_file():
        raise BenchmarkError(f"Missing timer report: {path}")
    sections: dict[str, dict[str, Any]] = {"inclusive": {}, "exclusive": {}}
    current: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped == "INCLUSIVE TIME:":
            current = "inclusive"
            continue
        if stripped == "EXCLUSIVE TIME:":
            current = "exclusive"
            continue
        if current is None or not re.match(r"^\s*\d+", line):
            continue
        try:
            timer_id = int(line[:6])
        except ValueError:
            continue
        if timer_id not in timer_ids:
            continue
        region = line[7:39].strip()
        values = line[40:].split()
        if len(values) != 5:
            raise BenchmarkError(f"Cannot parse timer line: {line}")
        sections[current][str(timer_id)] = {
            "region": region,
            "average": float(values[0]),
            "minimum": float(values[1]),
            "maximum": float(values[2]),
            "instances": float(values[3]),
            "total": float(values[4]),
        }
    if not sections["inclusive"]:
        raise BenchmarkError(f"No requested timers found in {path}")
    return sections


def require_hdf5_support() -> None:
    if h5py is None or np is None:
        raise BenchmarkError(
            "Correctness checks require Python packages h5py and numpy"
        )


def read_particles(path: Path) -> dict[str, Any]:
    require_hdf5_support()
    if not path.is_file():
        raise BenchmarkError(f"Missing particle file: {path}")
    with h5py.File(path, "r") as handle:
        missing = [name for name in ("ipart", "x", "y", "z") if name not in handle]
        if missing:
            raise BenchmarkError(f"{path} lacks datasets: {', '.join(missing)}")
        ids = np.asarray(handle["ipart"]).reshape(-1)
        coordinates = np.column_stack(
            [np.asarray(handle[name]).reshape(-1) for name in ("x", "y", "z")]
        )
        optional_fields = {}
        for name in ("state", "igrid", "seed"):
            if name in handle:
                optional_fields[name] = np.asarray(handle[name]).reshape(-1)
        schema = handle.attrs.get("PARTICLE_SCHEMA_VERSION")
    if len(ids) != len(coordinates):
        raise BenchmarkError(f"Particle dataset lengths differ in {path}")
    if len(np.unique(ids)) != len(ids):
        raise BenchmarkError(f"Duplicate particle IDs in {path}")
    if not np.isfinite(coordinates).all():
        raise BenchmarkError(f"Non-finite particle coordinates in {path}")
    order = np.argsort(ids, kind="stable")
    schema_value = int(schema) if schema is not None else None
    if schema_value is not None:
        for required in ("state", "igrid"):
            if required not in optional_fields:
                raise BenchmarkError(
                    f"Schema-versioned particle file {path} lacks {required}"
                )
    for name, values in optional_fields.items():
        if len(values) != len(ids):
            raise BenchmarkError(f"Particle dataset {name} has wrong length in {path}")
    return {
        "ids": ids[order],
        "coordinates": coordinates[order],
        "schema_version": schema_value,
        "fields": {name: values[order] for name, values in optional_fields.items()},
    }


def validate_particle_file(
    path: Path,
    expected_particles: int,
    expected_coordinates: list[list[float]] | None = None,
    tolerance: float = 0.0,
    coordinate_bounds: dict[str, list[float]] | None = None,
) -> dict[str, Any]:
    particles = read_particles(path)
    count = len(particles["ids"])
    if count != expected_particles:
        raise BenchmarkError(
            f"Expected {expected_particles} particles in {path}, found {count}"
        )
    result = {
        "particle_count": count,
        "minimum_id": int(particles["ids"][0]) if count else None,
        "maximum_id": int(particles["ids"][-1]) if count else None,
        "schema_version": particles["schema_version"],
    }
    if expected_coordinates is not None:
        expected = np.asarray(expected_coordinates, dtype=float)
        if expected.shape != particles["coordinates"].shape:
            raise BenchmarkError(f"Expected coordinate shape is wrong for {path}")
        max_delta = float(
            np.abs(particles["coordinates"] - expected).max(initial=0.0)
        )
        if max_delta > tolerance:
            raise BenchmarkError(
                f"Coordinate delta {max_delta:.6e} exceeds {tolerance:.6e} in {path}"
            )
        result["maximum_coordinate_delta"] = max_delta
        result["coordinate_tolerance"] = tolerance
    if coordinate_bounds is not None:
        minimum = np.asarray(coordinate_bounds["minimum"], dtype=float)
        maximum = np.asarray(coordinate_bounds["maximum"], dtype=float)
        if minimum.shape != (3,) or maximum.shape != (3,) or np.any(minimum > maximum):
            raise BenchmarkError(f"Coordinate bounds are invalid for {path}")
        coordinates = particles["coordinates"]
        if np.any(coordinates < minimum) or np.any(coordinates > maximum):
            raise BenchmarkError(
                f"Particle coordinates fall outside configured bounds in {path}"
            )
        result["coordinate_bounds"] = {
            "minimum": minimum.tolist(),
            "maximum": maximum.tolist(),
        }
    return result


def compare_particle_files(
    reference: Path,
    candidate: Path,
    expected_particles: int,
    tolerance: float,
) -> dict[str, Any]:
    reference_particles = read_particles(reference)
    candidate_particles = read_particles(candidate)
    if len(reference_particles["ids"]) != expected_particles:
        raise BenchmarkError("Reference particle count is incorrect")
    if len(candidate_particles["ids"]) != expected_particles:
        raise BenchmarkError("Restarted particle count is incorrect")
    if not np.array_equal(reference_particles["ids"], candidate_particles["ids"]):
        raise BenchmarkError("Particle IDs differ after restart")
    if reference_particles["schema_version"] != candidate_particles["schema_version"]:
        raise BenchmarkError("Particle schema version differs after restart")
    if reference_particles["fields"].keys() != candidate_particles["fields"].keys():
        raise BenchmarkError("Stored particle fields differ after restart")
    for name in reference_particles["fields"]:
        if not np.array_equal(
            reference_particles["fields"][name],
            candidate_particles["fields"][name],
        ):
            raise BenchmarkError(f"Stored particle field {name} differs after restart")
    coordinate_delta = np.abs(
        reference_particles["coordinates"] - candidate_particles["coordinates"]
    )
    max_delta = float(coordinate_delta.max(initial=0.0))
    if max_delta > tolerance:
        raise BenchmarkError(
            f"Restart coordinate delta {max_delta:.6e} exceeds {tolerance:.6e}"
        )
    return {
        "particle_count": expected_particles,
        "maximum_coordinate_delta": max_delta,
        "coordinate_tolerance": tolerance,
    }


def run_smoke(
    entry: dict[str, Any],
    root: Path,
    binary: Path,
    mpirun: str,
    mpirun_args: list[str],
) -> dict[str, Any]:
    run_dir = root / entry["name"]
    stage_case(entry, run_dir)
    run = run_mglet(binary, run_dir, entry["ranks"], mpirun, mpirun_args)
    validation = validate_particle_file(
        run_dir / "particles.h5",
        entry["expected_particles"],
        entry.get("expected_coordinates"),
        float(entry.get("coordinate_tolerance", 0.0)),
        entry.get("coordinate_bounds"),
    )
    return {"name": entry["name"], "kind": "smoke", "run": run, **validation}


def restart_entry(
    entry: dict[str, Any],
    steps: int,
    read_h5: bool,
) -> dict[str, Any]:
    variant = copy.deepcopy(entry)
    deep_merge(
        variant.setdefault("overrides", {}),
        {
            "time": {"mtstep": steps},
            "particles": {
                "dread_part_h5": read_h5,
                "dread_part_dict": not read_h5,
                "dwrite_part_h5": True,
            },
        },
    )
    return variant


def run_restart(
    entry: dict[str, Any],
    root: Path,
    binary: Path,
    mpirun: str,
    mpirun_args: list[str],
) -> dict[str, Any]:
    total_steps = int(entry["total_steps"])
    split_steps = int(entry["split_steps"])
    remaining_steps = total_steps - split_steps
    if split_steps <= 0 or remaining_steps <= 0:
        raise BenchmarkError(f"Invalid restart step split for {entry['name']}")

    case_root = root / entry["name"]
    reference_dir = case_root / "reference"
    first_dir = case_root / "segment-a"
    second_dir = case_root / "segment-b"

    stage_case(restart_entry(entry, total_steps, False), reference_dir)
    reference_run = run_mglet(
        binary, reference_dir, entry["ranks"], mpirun, mpirun_args
    )
    stage_case(restart_entry(entry, split_steps, False), first_dir)
    first_run = run_mglet(binary, first_dir, entry["ranks"], mpirun, mpirun_args)
    stage_case(restart_entry(entry, remaining_steps, True), second_dir)
    shutil.copy2(first_dir / "particles.h5", second_dir / "particles.h5")
    second_run = run_mglet(binary, second_dir, entry["ranks"], mpirun, mpirun_args)

    comparison = compare_particle_files(
        reference_dir / "particles.h5",
        second_dir / "particles.h5",
        entry["expected_particles"],
        float(entry["coordinate_tolerance"]),
    )
    return {
        "name": entry["name"],
        "kind": "restart",
        "reference_run": reference_run,
        "segment_a_run": first_run,
        "segment_b_run": second_run,
        **comparison,
    }


def run_correctness(
    entries: list[dict[str, Any]],
    root: Path,
    binary: Path,
    mpirun: str,
    mpirun_args: list[str],
) -> list[dict[str, Any]]:
    results = []
    for entry in entries:
        print(f"[correctness] {entry['name']}", flush=True)
        if entry["kind"] == "smoke":
            result = run_smoke(entry, root, binary, mpirun, mpirun_args)
        elif entry["kind"] == "restart":
            result = run_restart(entry, root, binary, mpirun, mpirun_args)
        else:
            raise BenchmarkError(f"Unknown correctness kind: {entry['kind']}")
        results.append(result)
    return results


def median_timer_runs(runs: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for section in ("inclusive", "exclusive"):
        timer_ids = sorted(
            {
                timer_id
                for run in runs
                for timer_id in run["timers"][section].keys()
            },
            key=int,
        )
        section_result: dict[str, Any] = {}
        for timer_id in timer_ids:
            available = [
                run["timers"][section][timer_id]
                for run in runs
                if timer_id in run["timers"][section]
            ]
            section_result[timer_id] = {
                "region": available[0]["region"],
                **{
                    field: statistics.median(item[field] for item in available)
                    for field in (
                        "average",
                        "minimum",
                        "maximum",
                        "instances",
                        "total",
                    )
                },
            }
        result[section] = section_result
    return result


def validate_timer_consistency(
    runs: list[dict[str, Any]], required_timer_ids: set[str]
) -> None:
    if not runs:
        raise BenchmarkError("Benchmark produced no measured runs")
    reference = {
        section: {
            timer_id: timer["region"]
            for timer_id, timer in runs[0]["timers"][section].items()
        }
        for section in ("inclusive", "exclusive")
    }
    for section in ("inclusive", "exclusive"):
        missing = required_timer_ids - reference[section].keys()
        if missing:
            raise BenchmarkError(
                f"Mandatory {section} timers missing: {', '.join(sorted(missing))}"
            )
    for index, run in enumerate(runs[1:], start=2):
        current = {
            section: {
                timer_id: timer["region"]
                for timer_id, timer in run["timers"][section].items()
            }
            for section in ("inclusive", "exclusive")
        }
        if current != reference:
            raise BenchmarkError(
                f"Timer IDs or region names changed in measured run {index}"
            )


def parse_initialized_particles(log_path: Path) -> int:
    pattern = re.compile(r"INITIALIZATION OF\s+(\d+)\s+PARTICLE\(S\)")
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            return int(match.group(1))
    raise BenchmarkError(f"Could not find initialized particle count in {log_path}")


def validate_benchmark_initialization(
    entry: dict[str, Any],
    case_root: Path,
    binary: Path,
    mpirun: str,
    mpirun_args: list[str],
) -> dict[str, Any]:
    preflight = copy.deepcopy(entry)
    deep_merge(
        preflight.setdefault("overrides", {}),
        {
            "time": {"mtstep": 0},
            "particles": {
                "terminal": "normal",
                "dwrite_part_h5": False,
                "snapshot_step": -1,
            },
        },
    )
    run_dir = case_root / "preflight"
    stage_case(preflight, run_dir)
    run = run_mglet(binary, run_dir, entry["ranks"], mpirun, mpirun_args)
    actual = parse_initialized_particles(run_dir / "mglet.OUT")
    expected = int(entry["configured_particles"])
    if actual != expected:
        raise BenchmarkError(
            f"{entry['name']} initialized {actual} particles, expected {expected}"
        )
    return {"expected_particles": expected, "initialized_particles": actual, "run": run}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_hashes(
    entry: dict[str, Any], parameters: dict[str, Any]
) -> dict[str, str]:
    source = source_path(entry)
    names = selected_h5_inputs(parameters)
    ib_config = parameters.get("ib", {})
    for geometry in ib_config.get("geometries", []):
        geometry_file = geometry.get("file")
        if isinstance(geometry_file, str):
            names.add(geometry_file)
    particles = parameters.get("particles", {})
    if particles.get("dread_part_dict", False):
        names.add("ParticleDict.txt")
    if particles.get("dread_obst_dict", False):
        names.add("ObstaclesDict.txt")
    external_paths = external_input_paths(entry)
    hashes = {}
    for name in sorted(names):
        path = external_paths.get(name, source / name)
        if not path.is_file():
            raise BenchmarkError(f"Required workload input is missing: {path}")
        hashes[name] = file_sha256(path)
    return hashes


def workload_fingerprint(
    entry: dict[str, Any], parameters: dict[str, Any]
) -> tuple[str, dict[str, Any]]:
    workload = {
        "source": source_identifier(entry),
        "staging": {
            "particle_dict_columns": entry.get(
                "particle_dict_columns", "id_xyz"
            )
        },
        "ranks": entry["ranks"],
        "configured_particles": entry["configured_particles"],
        "parameters": parameters,
        "input_sha256": input_hashes(entry, parameters),
    }
    encoded = json.dumps(workload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest(), workload


def run_benchmark(
    entry: dict[str, Any],
    root: Path,
    binary: Path,
    timer_ids: set[int],
    mpirun: str,
    mpirun_args: list[str],
) -> dict[str, Any]:
    case_root = root / entry["name"]
    preflight = validate_benchmark_initialization(
        entry, case_root, binary, mpirun, mpirun_args
    )
    measured_runs: list[dict[str, Any]] = []
    effective_parameters: dict[str, Any] | None = None
    total_runs = int(entry["warmups"]) + int(entry["repetitions"])
    for index in range(total_runs):
        is_warmup = index < int(entry["warmups"])
        category = "warmup" if is_warmup else "run"
        category_index = index + 1 if is_warmup else index - int(entry["warmups"]) + 1
        run_dir = case_root / f"{category}-{category_index:02d}"
        print(f"[benchmark] {entry['name']} {category} {category_index}", flush=True)
        parameters = stage_case(entry, run_dir)
        effective_parameters = parameters
        run = run_mglet(binary, run_dir, entry["ranks"], mpirun, mpirun_args)
        run["timers"] = parse_timer_report(
            run_dir / "mglet-perf-report.txt", timer_ids
        )
        run["parameters"] = parameters
        run["snapshot_bytes"] = directory_size(run_dir / "Particle_Snapshots")
        if not is_warmup:
            measured_runs.append(run)
    required_timers = {
        str(timer_id) for timer_id in entry.get("required_timers", [920, 940])
    }
    validate_timer_consistency(measured_runs, required_timers)
    assert effective_parameters is not None
    fingerprint, workload = workload_fingerprint(entry, effective_parameters)
    particle_steps = int(entry["configured_particles"]) * int(
        effective_parameters["time"]["mtstep"]
    )
    return {
        "name": entry["name"],
        "tier": entry["tier"],
        "category": entry.get("category", "compute"),
        "ranks": entry["ranks"],
        "configured_particles": entry["configured_particles"],
        "particle_steps": particle_steps,
        "required_timers": sorted(required_timers, key=int),
        "warmups": entry["warmups"],
        "repetitions": entry["repetitions"],
        "preflight": preflight,
        "workload": workload,
        "workload_fingerprint": fingerprint,
        "runs": measured_runs,
        "medians": {
            "wall_seconds": statistics.median(
                run["wall_seconds"] for run in measured_runs
            ),
            "output_bytes": statistics.median(
                run["output_bytes"] for run in measured_runs
            ),
            "snapshot_bytes": statistics.median(
                run["snapshot_bytes"] for run in measured_runs
            ),
            "timers": median_timer_runs(measured_runs),
        },
    }


def cmake_metadata(binary: Path) -> dict[str, str]:
    for directory in (binary.parent, *binary.parents):
        cache = directory / "CMakeCache.txt"
        if not cache.is_file():
            continue
        wanted = {
            "CMAKE_BUILD_TYPE",
            "CMAKE_C_COMPILER",
            "CMAKE_CXX_COMPILER",
            "CMAKE_Fortran_COMPILER",
            "MGLET_C_FLAGS",
            "MGLET_CXX_FLAGS",
            "MGLET_Fortran_FLAGS",
            "MGLET_OPENMP",
            "MGLET_REAL64",
        }
        values: dict[str, str] = {}
        for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith(("#", "//")) or "=" not in line:
                continue
            key_and_type, value = line.split("=", 1)
            key = key_and_type.split(":", 1)[0]
            if key in wanted:
                values[key] = value
        values["cache"] = str(cache)
        return values
    return {"cache": "not found"}


def cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name") and ":" in line:
                return line.split(":", 1)[1].strip()
    return platform.processor() or "unavailable"


def git_metadata() -> dict[str, Any]:
    status = command_output(["git", "status", "--porcelain"], REPO_ROOT)
    return {
        "repository": str(REPO_ROOT),
        "branch": command_output(["git", "branch", "--show-current"], REPO_ROOT),
        "commit": command_output(["git", "rev-parse", "HEAD"], REPO_ROOT),
        "dirty": bool(status and status != "unavailable"),
        "status": status.splitlines() if status != "unavailable" else [],
    }


def environment_metadata(binary: Path, mpirun: str, mpirun_args: list[str]) -> dict[str, Any]:
    return {
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "cpu_model": cpu_model(),
        "logical_cpus": os.cpu_count(),
        "python": sys.version.split()[0],
        "binary": str(binary),
        "binary_size": binary.stat().st_size,
        "mpi_version": command_output([mpirun, "--version"]).splitlines()[0],
        "mpirun": mpirun,
        "mpirun_args": mpirun_args,
        "loaded_modules": os.environ.get("LOADEDMODULES", ""),
        "affinity_environment": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "OMP_PROC_BIND",
                "OMP_PLACES",
                "SLURM_CPU_BIND",
            )
        },
        "cmake": cmake_metadata(binary),
        "git": git_metadata(),
    }


def summary_lines(summary: dict[str, Any]) -> list[str]:
    lines = [
        "Particle benchmark summary",
        f"Session: {summary['session_id']}",
        f"Commit: {summary['environment']['git']['commit']}",
        f"Host: {summary['environment']['hostname']}",
        "",
        "Correctness:",
    ]
    if summary["correctness"]:
        for result in summary["correctness"]:
            lines.append(f"  PASS {result['name']}")
    else:
        lines.append("  SKIPPED")
    lines.extend(("", "Benchmarks (median):"))
    for result in summary["benchmarks"]:
        timer_920 = (
            result["medians"]["timers"]["inclusive"]
            .get("920", {})
            .get("total")
        )
        timer_940 = (
            result["medians"]["timers"]["inclusive"]
            .get("940", {})
            .get("total")
        )
        lines.append(
            f"  {result['name']} [{result['category']}]: "
            f"wall={result['medians']['wall_seconds']:.6g}s "
            f"timer920={timer_920!s} timer940={timer_940!s} "
            f"snapshots={result['medians']['snapshot_bytes']}B"
        )
    return lines


def select_benchmarks(
    entries: list[dict[str, Any]], tier: str, names: list[str]
) -> list[dict[str, Any]]:
    selected = [
        entry
        for entry in entries
        if (tier == "all" or entry["tier"] == tier)
        and (not names or entry["name"] in names)
    ]
    unknown = set(names) - {entry["name"] for entry in entries}
    if unknown:
        raise BenchmarkError(f"Unknown benchmark case(s): {', '.join(sorted(unknown))}")
    if not selected:
        raise BenchmarkError("No benchmark cases selected")
    return selected


def run_command(args: argparse.Namespace) -> int:
    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = load_manifest(manifest_path)
    binary = require_within_workspace(Path(args.binary), "MGLET binary")
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise BenchmarkError(f"MGLET binary is not executable: {binary}")
    results_root = require_within_workspace(Path(args.results_root), "Results root")
    session_id = dt.datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    session_root = results_root / session_id
    session_root.mkdir(parents=True, exist_ok=False)

    summary: dict[str, Any] = {
        "schema_version": 1,
        "session_id": session_id,
        "manifest": str(manifest_path),
        "pilot": bool(args.pilot),
        "environment": environment_metadata(binary, args.mpirun, args.mpirun_arg),
        "correctness": [],
        "benchmarks": [],
    }
    write_json(session_root / "environment.json", summary["environment"])
    try:
        if not args.skip_correctness:
            summary["correctness"] = run_correctness(
                manifest["correctness"],
                session_root / "correctness",
                binary,
                args.mpirun,
                args.mpirun_arg,
            )
        if not args.correctness_only:
            entries = select_benchmarks(manifest["benchmarks"], args.tier, args.case)
            for entry in entries:
                entry = copy.deepcopy(entry)
                if args.pilot:
                    entry["warmups"] = 0
                    entry["repetitions"] = 1
                summary["benchmarks"].append(
                    run_benchmark(
                        entry,
                        session_root / "benchmarks",
                        binary,
                        set(manifest["timers"]),
                        args.mpirun,
                        args.mpirun_arg,
                    )
                )
    except Exception:
        summary["status"] = "failed"
        write_json(session_root / "summary.json", summary)
        raise
    summary["status"] = "passed"
    write_json(session_root / "summary.json", summary)
    (session_root / "summary.txt").write_text(
        "\n".join(summary_lines(summary)) + "\n", encoding="utf-8"
    )
    print(f"Results: {session_root}")
    return 0


def baseline_command(args: argparse.Namespace) -> int:
    result_path = require_within_workspace(Path(args.result), "Result file")
    output_path = require_within_workspace(Path(args.output), "Baseline file")
    result = read_json(result_path)
    if (
        result.get("schema_version") != 1
        or result.get("status") != "passed"
        or result.get("pilot", False)
    ):
        raise BenchmarkError("Only a passed benchmark result can become a baseline")
    baseline = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_result": str(result_path),
        "environment": result["environment"],
        "benchmarks": result["benchmarks"],
    }
    write_json(output_path, baseline)
    print(f"Baseline: {output_path}")
    return 0


def percentage_delta(current: float, baseline: float) -> float | None:
    if baseline == 0:
        return None
    return 100.0 * (current - baseline) / baseline


def comparable_environment_value(environment: dict[str, Any], key: str) -> Any:
    value = environment.get(key)
    if key == "cmake" and isinstance(value, dict):
        value = {item_key: item for item_key, item in value.items() if item_key != "cache"}
    return value


def compare_command(args: argparse.Namespace) -> int:
    result_path = require_within_workspace(Path(args.result), "Result file")
    baseline_path = require_within_workspace(Path(args.baseline), "Baseline file")
    result = read_json(result_path)
    baseline = read_json(baseline_path)
    if result.get("schema_version") != 1 or result.get("status") != "passed":
        raise BenchmarkError("Current result must be a passed schema-version 1 result")
    if baseline.get("schema_version") != 1:
        raise BenchmarkError("Baseline must use schema version 1")
    baseline_cases = {case["name"]: case for case in baseline["benchmarks"]}
    environment_keys = (
        "hostname",
        "platform",
        "cpu_model",
        "logical_cpus",
        "mpi_version",
        "mpirun_args",
        "loaded_modules",
        "affinity_environment",
        "cmake",
    )
    environment_match = {
        key: comparable_environment_value(result["environment"], key)
        == comparable_environment_value(baseline["environment"], key)
        for key in environment_keys
    }
    mismatches = [key for key, matches in environment_match.items() if not matches]
    if mismatches and not args.allow_environment_mismatch:
        raise BenchmarkError(
            "Environment mismatch; refusing comparison: "
            + ", ".join(mismatches)
            + ". Use --allow-environment-mismatch only for exploratory analysis."
        )
    comparisons = []
    for current in result["benchmarks"]:
        previous = baseline_cases.get(current["name"])
        if previous is None:
            continue
        if current.get("workload_fingerprint") != previous.get(
            "workload_fingerprint"
        ):
            raise BenchmarkError(
                f"Workload mismatch for {current['name']}; refusing comparison"
            )
        for section in ("inclusive", "exclusive"):
            current_regions = {
                timer_id: timer["region"]
                for timer_id, timer in current["medians"]["timers"][section].items()
            }
            previous_regions = {
                timer_id: timer["region"]
                for timer_id, timer in previous["medians"]["timers"][section].items()
            }
            for timer_id, region in current_regions.items():
                previous_region = previous_regions.get(timer_id)
                if previous_region is None:
                    raise BenchmarkError(
                        f"{section.capitalize()} timer {timer_id} ({region}) is "
                        f"absent from the baseline for {current['name']}"
                    )
                if previous_region != region:
                    raise BenchmarkError(
                        f"{section.capitalize()} timer {timer_id} changed identity "
                        f"for {current['name']}"
                    )
        metrics = {
            "wall_seconds": {
                "baseline": previous["medians"]["wall_seconds"],
                "current": current["medians"]["wall_seconds"],
            }
        }
        for section in ("inclusive", "exclusive"):
            previous_section = previous["medians"]["timers"][section]
            current_section = current["medians"]["timers"][section]
            for timer_id, old_timer in previous_section.items():
                timer = current_section.get(timer_id)
                # Disabled kernels may omit a timer entirely; treat that as zero work.
                current_total = 0.0 if timer is None else timer["total"]
                metrics[f"{section}.timer_{timer_id}.total"] = {
                    "baseline": old_timer["total"],
                    "current": current_total,
                }
        for values in metrics.values():
            values["delta_percent"] = percentage_delta(
                values["current"], values["baseline"]
            )
        comparisons.append({"name": current["name"], "metrics": metrics})
    if not comparisons:
        raise BenchmarkError("Result and baseline have no benchmark cases in common")
    comparison = {
        "schema_version": 1,
        "result": str(result_path),
        "baseline": str(baseline_path),
        "environment_match": environment_match,
        "environment_mismatch_allowed": bool(args.allow_environment_mismatch),
        "comparisons": comparisons,
    }
    output_path = (
        require_within_workspace(Path(args.output), "Comparison file")
        if args.output
        else result_path.with_name("comparison.json")
    )
    write_json(output_path, comparison)
    for case in comparisons:
        wall = case["metrics"]["wall_seconds"]
        print(
            f"{case['name']}: wall {wall['baseline']:.6g}s -> "
            f"{wall['current']:.6g}s ({wall['delta_percent']:+.2f}%)"
        )
    print(f"Comparison: {output_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run correctness and benchmarks")
    run_parser.add_argument("--binary", required=True, help="path to MGLET executable")
    run_parser.add_argument(
        "--tier",
        choices=("quick", "medium", "overnight", "all"),
        default="quick",
    )
    run_parser.add_argument(
        "--case", action="append", default=[], help="benchmark name; repeatable"
    )
    run_parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    run_parser.add_argument("--results-root", default=str(DEFAULT_RESULTS_ROOT))
    run_parser.add_argument("--mpirun", default="mpirun")
    run_parser.add_argument(
        "--mpirun-arg",
        action="append",
        default=[],
        help="extra mpirun argument; use --mpirun-arg=VALUE",
    )
    run_parser.add_argument("--skip-correctness", action="store_true")
    run_parser.add_argument(
        "--correctness-only",
        action="store_true",
        help="run correctness gates without timed benchmark cases",
    )
    run_parser.add_argument(
        "--pilot",
        action="store_true",
        help="run one measured repetition and mark result ineligible as a baseline",
    )
    run_parser.set_defaults(function=run_command)

    baseline_parser = subparsers.add_parser(
        "record-baseline", help="save a passed summary as a baseline"
    )
    baseline_parser.add_argument("--result", required=True)
    baseline_parser.add_argument("--output", required=True)
    baseline_parser.set_defaults(function=baseline_command)

    compare_parser = subparsers.add_parser(
        "compare", help="compare a passed summary with a baseline"
    )
    compare_parser.add_argument("--result", required=True)
    compare_parser.add_argument("--baseline", required=True)
    compare_parser.add_argument("--output")
    compare_parser.add_argument(
        "--allow-environment-mismatch",
        action="store_true",
        help="allow exploratory comparison across different environments",
    )
    compare_parser.set_defaults(function=compare_command)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.function(args)
    except BenchmarkError as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
