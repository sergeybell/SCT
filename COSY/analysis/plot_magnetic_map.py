#!/usr/bin/env python3
"""Plot Delta nu_s vs target chromaticity (Map_lat methodology for magnetic_2..8)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "COSY" / "analysis"))

from map_lat_lib import (  # noqa: E402
    DEFAULT_ROOT,
    MAGNETIC_LATTICES,
    PANEL_CHROM_LABELS,
    PANEL_COLORS,
    PANEL_CONSTRAINTS,
    PLANE_NAMES,
    STEPS_DEFAULT,
    build_lattice_models,
    compute_panel_curves,
    get_data,
    standard_chrom_ranges,
)


def _style_axes(ax) -> None:
    ax.axhline(0, color="black", lw=1.0, alpha=0.5, zorder=0)
    ax.axvline(0, color="black", lw=1.0, alpha=0.5, zorder=0)
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((-2, 2))
    ax.yaxis.set_major_formatter(fmt)
    ax.grid(True, linestyle="--", alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)


def plot_lattice(
    lat: str,
    root: Path,
    *,
    steps: int = STEPS_DEFAULT,
    require_fr3: bool = False,
    show: bool = False,
) -> bool:
    df0 = get_data(lat, "FR0", root)
    df3 = get_data(lat, "FR3", root)
    if df0 is None:
        print(f"Skip {lat}: no FR0 data")
        return False
    if require_fr3 and df3 is None:
        print(f"Skip {lat}: no FR3 data")
        return False

    for df, mode in ((df0, "FR0"),) + ((df3, "FR3"),) if df3 is not None else ():
        for col in ("chrom_x", "chrom_y", "eta_1"):
            if df[col].nunique() < 2:
                print(f"WARNING {lat} {mode}: {col} is constant — check sextupole wiring in .fox")

    models = build_lattice_models(df0, df3)
    ranges, curves = compute_panel_curves(models, standard_chrom_ranges(steps))

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=False)
    title = f"$\\Delta\\nu_s$ vs target chromaticity: {lat}"
    if df3 is None:
        title += " (FR0 only)"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    out_dir = root / lat / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    for i in range(3):
        ax = axes[i]
        color = PANEL_COLORS[i]
        if "WithFr_P0" in models:
            ax.plot(ranges[i], curves["WithFr_P0"][i], color=color, ls="-", lw=2.5, label="With Fringe P0")
            ax.plot(ranges[i], curves["WithFr_P1"][i], color=color, ls="--", lw=2.5, label="With Fringe P1")
        ax.plot(ranges[i], curves["NoFr_P0"][i], color="gray", ls=":", alpha=0.8, lw=2.0, label="No Fringe P0")
        ax.plot(ranges[i], curves["NoFr_P1"][i], color="gray", ls="-.", alpha=0.8, lw=2.0, label="No Fringe P1")

        ax.set_xlabel(PANEL_CHROM_LABELS[i], fontsize=12)
        ax.set_ylabel(r"$\Delta\nu_s$", fontsize=12)
        ax.set_title(f"Plane {PLANE_NAMES[i].upper()}\n({PANEL_CONSTRAINTS[i]})", fontweight="bold", fontsize=10)
        _style_axes(ax)
        ax.legend(fontsize=8, loc="best")

        t_fig, t_ax = plt.subplots(figsize=(6, 5))
        if "WithFr_P0" in models:
            t_ax.plot(ranges[i], curves["WithFr_P0"][i], color=color, ls="-", lw=2.5, label="With Fringe P0")
            t_ax.plot(ranges[i], curves["WithFr_P1"][i], color=color, ls="--", lw=2.5, label="With Fringe P1")
        t_ax.plot(ranges[i], curves["NoFr_P0"][i], color="gray", ls=":", alpha=0.6, label="No Fringe P0")
        t_ax.plot(ranges[i], curves["NoFr_P1"][i], color="gray", ls="-.", alpha=0.6, label="No Fringe P1")
        t_ax.set_xlabel(PANEL_CHROM_LABELS[i])
        t_ax.set_ylabel(r"$\Delta\nu_s$")
        t_ax.set_title(f"{lat} | {PLANE_NAMES[i].upper()}-Plane | {PANEL_CONSTRAINTS[i]}", fontsize=11)
        _style_axes(t_ax)
        t_ax.legend(fontsize=9, loc="best")
        t_fig.savefig(out_dir / f"{lat}_compare_dual_control_{PLANE_NAMES[i]}.png", dpi=300, bbox_inches="tight")
        plt.close(t_fig)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    summary = out_dir / f"{lat}_dnu_s_vs_chrom.png"
    fig.savefig(summary, dpi=150, bbox_inches="tight", facecolor="white", edgecolor="none")
    if show:
        plt.show()
    else:
        plt.close(fig)
    print(f"OK: {summary}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot Delta nu_s vs chromaticity (Map_lat style)")
    ap.add_argument("lattices", nargs="*", default=MAGNETIC_LATTICES)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--require-fr3", action="store_true")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    ok = False
    for lat in args.lattices:
        if plot_lattice(lat, args.root, require_fr3=args.require_fr3, show=args.show):
            ok = True
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
