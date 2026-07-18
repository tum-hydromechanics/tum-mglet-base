#!/usr/bin/env python3
"""Unit tests for the particle benchmark harness."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SCRIPT = Path(__file__).with_name("particle_bench.py")
SPEC = importlib.util.spec_from_file_location("particle_bench", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
particle_bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(particle_bench)


class ParticleBenchTests(unittest.TestCase):
    def temporary_directory(self) -> tempfile.TemporaryDirectory[str]:
        return tempfile.TemporaryDirectory(dir=particle_bench.WORKSPACE_ROOT)

    def test_deep_merge_preserves_unmodified_values(self) -> None:
        target = {"time": {"mtstep": 100, "dt": 0.1}, "flow": {"solve": True}}
        particle_bench.deep_merge(
            target, {"time": {"mtstep": 20}, "particles": {"terminal": "none"}}
        )
        self.assertEqual(target["time"], {"mtstep": 20, "dt": 0.1})
        self.assertTrue(target["flow"]["solve"])
        self.assertEqual(target["particles"]["terminal"], "none")

    def test_parse_timer_report_reads_both_sections(self) -> None:
        def timer_line(timer_id: int, region: str, total: float) -> str:
            return (
                f"{timer_id:6d} {region:<32} "
                f"{1.0:12.3E}{0.9:12.3E}{1.1:12.3E}"
                f"{2.0:12.3E}{total:12.3E}"
            )

        report = "\n".join(
            [
                "INCLUSIVE TIME:",
                "   idx region",
                timer_line(920, "PSIM_TIMEINTEGRATION", 2.5),
                timer_line(940, "PSIM_EXCHANGE", 0.5),
                "",
                "EXCLUSIVE TIME:",
                "   idx region",
                timer_line(920, "PSIM_TIMEINTEGRATION", 2.0),
            ]
        )
        with self.temporary_directory() as directory:
            path = Path(directory) / "mglet-perf-report.txt"
            path.write_text(report, encoding="utf-8")
            parsed = particle_bench.parse_timer_report(path, {920, 940})
        self.assertEqual(parsed["inclusive"]["920"]["region"], "PSIM_TIMEINTEGRATION")
        self.assertEqual(parsed["inclusive"]["940"]["total"], 0.5)
        self.assertEqual(parsed["exclusive"]["920"]["total"], 2.0)

    def test_particle_comparison_sorts_ids(self) -> None:
        with self.temporary_directory() as directory:
            root = Path(directory)
            reference = root / "reference.h5"
            candidate = root / "candidate.h5"
            self.write_particles(reference, [2, 1], [0.2, 0.1])
            self.write_particles(candidate, [1, 2], [0.1, 0.2])
            comparison = particle_bench.compare_particle_files(
                reference, candidate, expected_particles=2, tolerance=0.0
            )
        self.assertEqual(comparison["particle_count"], 2)
        self.assertEqual(comparison["maximum_coordinate_delta"], 0.0)

    def test_safe_name_rejects_path_traversal(self) -> None:
        with self.assertRaises(particle_bench.BenchmarkError):
            particle_bench.safe_name("../outside", "case name")

    def test_parameter_paths_reject_absolute_output(self) -> None:
        with self.assertRaises(particle_bench.BenchmarkError):
            particle_bench.validate_parameter_paths(
                {"io": {"grids": "grid.h5", "outfile": "/tmp/fields.h5"}}
            )

    def test_manifests_and_external_bindings_are_valid(self) -> None:
        particle_bench.load_manifest(SCRIPT.with_name("benchmarks.json"))
        manifest = particle_bench.load_manifest(
            SCRIPT.with_name("bcc-benchmarks.json")
        )
        entry = manifest["benchmarks"][0]
        self.assertTrue(particle_bench.source_path(entry).is_dir())
        bindings = particle_bench.external_input_paths(entry)
        self.assertEqual(len(bindings), 4)
        self.assertTrue(all(path.is_file() for path in bindings.values()))

    def test_generated_seed_covers_every_rank(self) -> None:
        parameters = {"particles": {}}
        entry = {
            "ranks": 4,
            "seed_base": 100,
            "seed_words_per_rank": 8,
        }
        particle_bench.apply_generated_seed(entry, parameters)
        self.assertEqual(len(parameters["particles"]["particle_seed"]), 32)
        self.assertEqual(parameters["particles"]["particle_seed"][0], 100)

    def test_legacy_particle_dictionary_is_converted_during_staging(self) -> None:
        with self.temporary_directory() as directory:
            source = Path(directory) / "ParticleDict.txt"
            destination = Path(directory) / "converted.txt"
            source.write_text("2\n0.1 0.2 0.3\n0.4 0.5 0.6\n", encoding="utf-8")
            particle_bench.copy_case_input(
                source, destination, {"particle_dict_columns": "xyz"}
            )
            self.assertEqual(
                destination.read_text(encoding="utf-8"),
                "2\n1 0.1 0.2 0.3\n2 0.4 0.5 0.6\n",
            )

    def test_exchange_timer_is_required_when_requested(self) -> None:
        timer = {
            "region": "region",
            "average": 1.0,
            "minimum": 1.0,
            "maximum": 1.0,
            "instances": 1.0,
            "total": 1.0,
        }
        run = {
            "timers": {
                "inclusive": {"920": timer, "940": timer},
                "exclusive": {"920": timer, "940": timer},
            }
        }
        with self.assertRaises(particle_bench.BenchmarkError):
            particle_bench.validate_timer_consistency([run], {"920", "940", "942"})

    @staticmethod
    def write_particles(path: Path, ids: list[int], x_coordinates: list[float]) -> None:
        with h5py.File(path, "w") as handle:
            handle.attrs["PARTICLE_SCHEMA_VERSION"] = 1
            handle.create_dataset("ipart", data=np.asarray(ids, dtype=np.int32))
            handle.create_dataset("state", data=np.ones(len(ids), dtype=np.int32))
            handle.create_dataset("igrid", data=np.ones(len(ids), dtype=np.int32))
            handle.create_dataset("x", data=np.asarray(x_coordinates))
            handle.create_dataset("y", data=np.zeros(len(ids)))
            handle.create_dataset("z", data=np.zeros(len(ids)))


if __name__ == "__main__":
    unittest.main()
