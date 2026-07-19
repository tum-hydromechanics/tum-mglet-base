#!/usr/bin/env python3
"""Scientific plots of particle timer hotspots and shared-kernel plan targets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parents[3]
DEFAULT_SELF = (
    WORKSPACE_ROOT
    / "benchmark-results"
    / "20260719-145914-762784"
    / "summary.json"
)
DEFAULT_BCC = (
    WORKSPACE_ROOT
    / "benchmark-results"
    / "20260719-150533-923471"
    / "summary.json"
)
DEFAULT_OUT = WORKSPACE_ROOT / "benchmark-results" / "plots" / "shared-kernel-plan"

PLAN_TARGET_TIMERS = {
    "921": "ADV_VELOCITY",
    "922": "ADV_MOTION",
    "924": "DIF_RN_GENERATION",
    "925": "DIF_MOTION",
}
DEFERRED_EMPHASIS = {
    "941": "PREP_COMM",
    "942": "MPI_COMM",
    "944": "SORT_PLIST",
}

CASE_ORDER = [
    "diffusion-focused-medium",
    "advection-diffusion-medium",
    "advection-diffusion-exchange-medium",
    "bcc-compute-medium",
]


def load_cases(paths: list[Path]) -> dict[str, dict[str, Any]]:
    cases: dict[str, dict[str, Any]] = {}
    for path in paths:
        summary = json.loads(path.read_text(encoding="utf-8"))
        for entry in summary["benchmarks"]:
            cases[entry["name"]] = {
                "wall": float(entry["medians"]["wall_seconds"]),
                "exclusive": entry["medians"]["timers"]["exclusive"],
                "source": str(path),
            }
    return cases


def timer_rows(case: dict[str, Any]) -> list[tuple[str, str, float, float, str]]:
    wall = case["wall"]
    rows = []
    for timer_id, timer in case["exclusive"].items():
        total = float(timer["total"])
        share = 100.0 * total / wall if wall > 0 else 0.0
        region = str(timer["region"])
        if timer_id in PLAN_TARGET_TIMERS:
            role = "plan_target"
        elif timer_id in DEFERRED_EMPHASIS:
            role = "deferred"
        else:
            role = "other"
        rows.append((timer_id, region, total, share, role))
    rows.sort(key=lambda item: item[3], reverse=True)
    return rows


def short_case(name: str) -> str:
    return (
        name.replace("-medium", "")
        .replace("advection-diffusion", "adv-diff")
        .replace("diffusion-focused", "diffusion")
        .replace("bcc-compute", "bcc")
    )


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.pdf")
    fig.savefig(out_dir / f"{stem}.png", dpi=160)
    plt.close(fig)


def plot_exclusive_share(cases: dict[str, dict[str, Any]], out_dir: Path) -> None:
    selected = [name for name in CASE_ORDER if name in cases]
    # Top timers across all cases for a stable legend set.
    totals: dict[str, float] = {}
    labels: dict[str, str] = {}
    for name in selected:
        for timer_id, region, total, _share, _role in timer_rows(cases[name]):
            totals[timer_id] = totals.get(timer_id, 0.0) + total
            labels[timer_id] = region
    top_ids = [
        timer_id
        for timer_id, _ in sorted(totals.items(), key=lambda item: item[1], reverse=True)[:10]
    ]

    fig, axes = plt.subplots(len(selected), 1, figsize=(10, 2.2 * len(selected)), sharex=True)
    if len(selected) == 1:
        axes = [axes]
    colors = plt.cm.tab10(np.linspace(0, 1, len(top_ids)))
    color_map = {timer_id: colors[i] for i, timer_id in enumerate(top_ids)}

    for ax, name in zip(axes, selected):
        rows = {timer_id: share for timer_id, _r, _t, share, _role in timer_rows(cases[name])}
        values = [rows.get(timer_id, 0.0) for timer_id in top_ids]
        y = np.arange(len(top_ids))
        ax.barh(y, values, color=[color_map[t] for t in top_ids])
        ax.set_yticks(y)
        ax.set_yticklabels([labels[t] for t in top_ids], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Exclusive time / wall (%)")
        ax.set_title(f"{name}  (wall={cases[name]['wall']:.2f}s)")
        ax.set_xlim(0, max(55, max(values) * 1.15 if values else 1))
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("Exclusive timer share of wall time", y=1.01)
    save_figure(fig, out_dir, "exclusive_timer_share_by_case")


def plot_plan_targets_overlay(cases: dict[str, dict[str, Any]], out_dir: Path) -> None:
    selected = [name for name in CASE_ORDER if name in cases]
    role_colors = {"plan_target": "#1f77b4", "deferred": "#ff7f0e", "other": "#7f7f7f"}
    fig, axes = plt.subplots(len(selected), 1, figsize=(10, 2.4 * len(selected)), sharex=True)
    if len(selected) == 1:
        axes = [axes]

    for ax, name in zip(axes, selected):
        rows = timer_rows(cases[name])[:12]
        shares = [share for _id, _r, _t, share, _role in rows]
        labels = [region for _id, region, _t, _s, _role in rows]
        colors = [role_colors[role] for _id, _r, _t, _s, role in rows]
        y = np.arange(len(rows))
        ax.barh(y, shares, color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(f"{name}  (wall={cases[name]['wall']:.2f}s)")
        ax.set_xlabel("Exclusive time / wall (%)")
        ax.grid(axis="x", alpha=0.3)
        ax.set_xlim(0, max(55, max(shares) * 1.15 if shares else 1))

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=role_colors["plan_target"], label="In this plan (921/922/924/925)"),
        plt.Rectangle((0, 0), 1, 1, color=role_colors["deferred"], label="Deferred (exchange/sort)"),
        plt.Rectangle((0, 0), 1, 1, color=role_colors["other"], label="Other"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Plan targets vs deferred exclusive costs", y=1.06)
    save_figure(fig, out_dir, "plan_targets_overlay")


def plot_absolute_targets(cases: dict[str, dict[str, Any]], out_dir: Path) -> None:
    selected = [name for name in CASE_ORDER if name in cases]
    timer_ids = list(PLAN_TARGET_TIMERS.keys())
    x = np.arange(len(selected))
    width = 0.18
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, timer_id in enumerate(timer_ids):
        values = []
        for name in selected:
            timer = cases[name]["exclusive"].get(timer_id)
            values.append(float(timer["total"]) if timer else 0.0)
        ax.bar(x + (i - 1.5) * width, values, width, label=f"{timer_id} {PLAN_TARGET_TIMERS[timer_id]}")
    ax.set_xticks(x)
    ax.set_xticklabels([short_case(name) for name in selected])
    ax.set_ylabel("Exclusive seconds")
    ax.set_title("Absolute exclusive time for plan-target timers")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    save_figure(fig, out_dir, "absolute_exclusive_seconds")


def plot_wall_breakdown_stacked(cases: dict[str, dict[str, Any]], out_dir: Path) -> None:
    selected = [name for name in CASE_ORDER if name in cases]
    categories = ["plan_targets", "deferred_exchange", "other_exclusive", "residual"]
    colors = {
        "plan_targets": "#1f77b4",
        "deferred_exchange": "#ff7f0e",
        "other_exclusive": "#7f7f7f",
        "residual": "#c7c7c7",
    }
    stacks = {key: [] for key in categories}
    for name in selected:
        wall = cases[name]["wall"]
        plan = 0.0
        deferred = 0.0
        other = 0.0
        for timer_id, timer in cases[name]["exclusive"].items():
            total = float(timer["total"])
            if timer_id in PLAN_TARGET_TIMERS:
                plan += total
            elif timer_id in DEFERRED_EMPHASIS:
                deferred += total
            else:
                other += total
        accounted = plan + deferred + other
        residual = max(0.0, wall - accounted)
        stacks["plan_targets"].append(plan)
        stacks["deferred_exchange"].append(deferred)
        stacks["other_exclusive"].append(other)
        stacks["residual"].append(residual)

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(selected))
    bottom = np.zeros(len(selected))
    labels = {
        "plan_targets": "Plan targets (921/922/924/925)",
        "deferred_exchange": "Deferred exchange/sort (941/942/944)",
        "other_exclusive": "Other exclusive timers",
        "residual": "Residual (wall − exclusive sum)",
    }
    for key in categories:
        values = np.asarray(stacks[key])
        ax.bar(x, values, bottom=bottom, color=colors[key], label=labels[key])
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels([short_case(name) for name in selected])
    ax.set_ylabel("Seconds")
    ax.set_title("Wall-time breakdown with plan-target segment")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    save_figure(fig, out_dir, "wall_breakdown_stacked")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-summary", type=Path, default=DEFAULT_SELF)
    parser.add_argument("--bcc-summary", type=Path, default=DEFAULT_BCC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    paths = [args.self_summary, args.bcc_summary]
    for path in paths:
        if not path.is_file():
            raise SystemExit(f"Missing summary: {path}")

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cases = load_cases(paths)

    plot_exclusive_share(cases, out_dir)
    plot_plan_targets_overlay(cases, out_dir)
    plot_absolute_targets(cases, out_dir)
    plot_wall_breakdown_stacked(cases, out_dir)

    manifest = {
        "inputs": [str(path) for path in paths],
        "cases": sorted(cases),
        "plan_target_timers": PLAN_TARGET_TIMERS,
        "deferred_timers": DEFERRED_EMPHASIS,
        "outputs": [
            "exclusive_timer_share_by_case.pdf",
            "exclusive_timer_share_by_case.png",
            "plan_targets_overlay.pdf",
            "plan_targets_overlay.png",
            "absolute_exclusive_seconds.pdf",
            "absolute_exclusive_seconds.png",
            "wall_breakdown_stacked.pdf",
            "wall_breakdown_stacked.png",
        ],
    }
    (out_dir / "plot_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Plots written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
