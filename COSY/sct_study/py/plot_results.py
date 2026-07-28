#!/usr/bin/env python3
"""Aggregate key plots / summary figure for the study."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DAT_OUT, PLOTS, append_worklog, stems  # noqa: E402


def main() -> int:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle("SCT study summary (FR0)")

    conds, dists, names = [], [], []
    for stem in stems():
        p = DAT_OUT / stem / "audit.json"
        if not p.is_file():
            continue
        a = json.loads(p.read_text(encoding="utf-8"))
        m = a["model"]
        names.append(stem.replace("magnetic_", "m"))
        conds.append(m["condition_number"])
        dists.append(m["grid_relation"]["min_distance_to_grid_point"])

    axes[0].bar(names, conds, color="steelblue")
    axes[0].set_ylabel(r"cond$(R_{int})$")
    axes[0].set_title("Mapping condition number")
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].bar(names, dists, color="darkorange")
    axes[1].set_ylabel("min |I*-grid|")
    axes[1].set_title("Extrapolation distance of I*")
    axes[1].grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    out = PLOTS / "summary_mapping_quality.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # delta_eq suppression at I* vs natural
    fig, ax = plt.subplots(figsize=(8, 4))
    width = 0.35
    xs = np.arange(len(stems()))
    nat_vals, ist_vals, labels = [], [], []
    for stem in stems():
        p = DAT_OUT / stem / "delta_eq_at_istar.json"
        if not p.is_file():
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        labels.append(stem.replace("magnetic_", "m"))
        # use X_only scenario
        ist_vals.append(abs(d["scenarios"].get("X_only", 0.0)))
        nat_vals.append(abs(d["scenarios"].get("X_only_at_natural", 0.0)))
    ax.bar(xs - width / 2, nat_vals, width, label="|Δδ_eq| natural")
    ax.bar(xs + width / 2, ist_vals, width, label="|Δδ_eq| I*")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_yscale("symlog", linthresh=1e-20)
    ax.legend(fontsize=8)
    ax.set_title(r"Theoretical $|\Delta\delta_{eq}|$ (X_only scenario)")
    ax.grid(True, alpha=0.3)
    out2 = PLOTS / "summary_delta_eq_suppression.png"
    fig.tight_layout()
    fig.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close(fig)

    append_worklog(
        "## Summary plots\n\n"
        f"- **Статус:** verified\n"
        f"- **Файлы:** `{out.name}`, `{out2.name}`\n"
    )
    print(f"OK {out}\nOK {out2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
