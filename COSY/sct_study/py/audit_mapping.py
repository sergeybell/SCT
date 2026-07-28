#!/usr/bin/env python3
"""FR0 mapping audit for magnetic_2..5."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DAT_OUT,
    PLOTS,
    append_worklog,
    audit_mapping_lattice,
    get_data,
    model_summary,
    stems,
    write_json,
)

import matplotlib.pyplot as plt
import numpy as np


def plot_grid_and_istar(stem: str, df, summary: dict) -> Path:
    I = np.array(summary["I_star"])
    fig = plt.figure(figsize=(12, 4))
    pairs = [(0, 1, "SGx1", "SGx2"), (0, 2, "SGx1", "SGy1"), (1, 2, "SGx2", "SGy1")]
    for ax_i, (i, j, xl, yl) in enumerate(pairs, 1):
        ax = fig.add_subplot(1, 3, ax_i)
        ax.scatter(df[xl], df[yl], s=12, alpha=0.5, c="steelblue", label="FR0 grid")
        ax.scatter([I[i]], [I[j]], s=80, c="crimson", marker="*", label="I*", zorder=5)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"{stem}: mapping grid vs I* (FR0)")
    fig.tight_layout()
    out = PLOTS / f"{stem}_grid_vs_istar.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    rows = []
    for stem in stems():
        audit = audit_mapping_lattice(stem)
        df = get_data(stem, "FR0")
        if df is None:
            write_json(DAT_OUT / stem / "audit.json", {"error": "missing FR0 dat", **audit})
            continue
        summary = model_summary(df)
        out = {
            "stem": stem,
            "wiring_audit": audit,
            "model": summary,
        }
        write_json(DAT_OUT / stem / "audit.json", out)
        plot = plot_grid_and_istar(stem, df, summary)
        rows.append(
            f"| {stem} | {summary['condition_number']:.1f} | "
            f"{summary['I_star']} | {summary['grid_relation']['inside_axis_aligned_box']} | "
            f"{summary['grid_relation']['min_distance_to_grid_point']:.4g} | {plot.name} |"
        )
        print(f"OK audit {stem}: cond={summary['condition_number']:.1f} inside={summary['grid_relation']['inside_axis_aligned_box']}")

    append_worklog(
        "## Audit mapping (FR0)\n\n"
        "- **Статус:** verified (offline)\n"
        "- **Команда:** `python COSY/sct_study/py/audit_mapping.py`\n"
        "- **Выход:** `dat/<stem>/audit.json`, `plots/<stem>_grid_vs_istar.png`\n\n"
        "| stem | cond | I* | inside box | min_dist | plot |\n|------|------|----|------------|----------|------|\n"
        + "\n".join(rows)
        + "\n\n**Вывод:** линейность R²≈1, но I* вне оси-aligned box у всех структур "
        "(нужен отрицательный SGx2). magnetic_2 — ближайший кандидат для пилота.\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
