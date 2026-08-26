#!/usr/bin/env python3
"""Mapping audit for FR0 / FR3 (presentation stems or study stems)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DAT_OUT,
    PLOTS,
    append_worklog,
    audit_mapping_lattice,
    get_data,
    mapping_presentation_stems,
    model_summary,
    stems,
    write_json,
)

import matplotlib.pyplot as plt
import numpy as np


def plot_grid_and_istar(stem: str, df, summary: dict, mode: str) -> Path:
    I = np.array(summary["I_star"])
    fig = plt.figure(figsize=(12, 4))
    pairs = [(0, 1, "SGx1", "SGx2"), (0, 2, "SGx1", "SGy1"), (1, 2, "SGx2", "SGy1")]
    for ax_i, (i, j, xl, yl) in enumerate(pairs, 1):
        ax = fig.add_subplot(1, 3, ax_i)
        ax.scatter(df[xl], df[yl], s=12, alpha=0.5, c="steelblue", label=f"{mode} grid")
        ax.scatter([I[i]], [I[j]], s=80, c="crimson", marker="*", label="I*", zorder=5)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"{stem}: mapping grid vs I* ({mode})")
    fig.tight_layout()
    out = PLOTS / f"{stem}_grid_vs_istar_{mode}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    # Keep legacy FR0 filename for magnetic study scripts
    if mode == "FR0":
        legacy = PLOTS / f"{stem}_grid_vs_istar.png"
        fig2 = plt.figure(figsize=(12, 4))
        for ax_i, (i, j, xl, yl) in enumerate(pairs, 1):
            ax = fig2.add_subplot(1, 3, ax_i)
            ax.scatter(df[xl], df[yl], s=12, alpha=0.5, c="steelblue", label="FR0 grid")
            ax.scatter([I[i]], [I[j]], s=80, c="crimson", marker="*", label="I*", zorder=5)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)
        fig2.suptitle(f"{stem}: mapping grid vs I* (FR0)")
        fig2.tight_layout()
        fig2.savefig(legacy, dpi=150, bbox_inches="tight")
        plt.close(fig2)
    return out


def audit_one(stem: str, mode: str) -> dict | None:
    audit = audit_mapping_lattice(stem)
    df = get_data(stem, mode)
    if df is None:
        write_json(
            DAT_OUT / stem / f"audit_{mode}.json",
            {"error": f"missing {mode} dat", **audit},
        )
        if mode == "FR0":
            write_json(DAT_OUT / stem / "audit.json", {"error": f"missing {mode} dat", **audit})
        return None
    summary = model_summary(df, stem=stem, mode=mode)
    out = {
        "stem": stem,
        "mode": mode,
        "wiring_audit": audit,
        "model": summary,
    }
    write_json(DAT_OUT / stem / f"audit_{mode}.json", out)
    if mode == "FR0":
        write_json(DAT_OUT / stem / "audit.json", out)
    plot = plot_grid_and_istar(stem, df, summary, mode)
    gr = summary["grid_relation"]
    d_norm = gr.get("min_distance_normalized", float("nan"))
    print(
        f"OK audit {stem} {mode}: cond={summary['condition_number']:.1f} "
        f"inside={gr['inside_axis_aligned_box']} d_norm={d_norm:.3g} plot={plot.name}"
    )
    return {
        "stem": stem,
        "mode": mode,
        "cond": summary["condition_number"],
        "I_star": summary["I_star"],
        "inside": gr["inside_axis_aligned_box"],
        "d_min": gr["min_distance_to_grid_point"],
        "d_norm": d_norm,
        "plot": plot.name,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Audit mapping grids (FR0/FR3)")
    ap.add_argument("--mode", choices=["FR0", "FR3", "both"], default="FR0")
    ap.add_argument("--stem", action="append", dest="stems", metavar="STEM")
    ap.add_argument(
        "--presentation",
        action="store_true",
        help="Use mapping_presentation_stems from study_config.json",
    )
    args = ap.parse_args()

    if args.stems:
        stem_list = list(args.stems)
    elif args.presentation:
        stem_list = mapping_presentation_stems()
    else:
        stem_list = stems()

    modes = ["FR0", "FR3"] if args.mode == "both" else [args.mode]
    rows = []
    for stem in stem_list:
        for mode in modes:
            row = audit_one(stem, mode)
            if row is None:
                continue
            rows.append(
                f"| {row['stem']} | {row['mode']} | {row['cond']:.1f} | "
                f"{row['I_star']} | {row['inside']} | "
                f"{row['d_min']:.4g} | {row['d_norm']:.3g} | {row['plot']} |"
            )

    append_worklog(
        "## Audit mapping\n\n"
        f"- **Статус:** verified (offline)\n"
        f"- **Команда:** `python COSY/sct_study/py/audit_mapping.py "
        f"--mode {args.mode}"
        + (" --presentation" if args.presentation else "")
        + "`\n"
        f"- **Выход:** `dat/<stem>/audit[_FR0|_FR3].json`, `plots/<stem>_grid_vs_istar*.png`\n\n"
        "| stem | mode | cond | I* | inside box | d_min | d_min/‖ΔI‖ | plot |\n"
        "|------|------|------|----|------------|-------|------------|------|\n"
        + "\n".join(rows)
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
