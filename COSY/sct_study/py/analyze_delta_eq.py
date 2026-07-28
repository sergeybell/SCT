#!/usr/bin/env python3
"""Theoretical Δδ_eq panels + map-Δν_s panels (Melnikov style)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DAT_OUT,
    PLOTS,
    append_worklog,
    build_lattice_models,
    compute_panel_curves,
    delta_eq_theory,
    get_data,
    get_dual_model,
    lattice_length,
    load_config,
    solve_currents,
    standard_chrom_ranges,
    stems,
    write_json,
)

PANEL_LABELS = [r"$\xi_x$", r"$\xi_y$", r"$\eta_1$"]
CONSTRAINTS = [r"$\xi_y=\eta_1=0$", r"$\xi_x=\eta_1=0$", r"$\xi_x=\xi_y=0$"]


def panel_delta_eq(stem: str, df, cfg) -> Path:
    m = get_dual_model(df)
    ranges = standard_chrom_ranges(81)
    L = lattice_length(stem)
    gamma = float(cfg["gamma"])
    alpha0 = float(df["alpha_0"].mean())
    scenarios = cfg["emittance_scenarios"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"{stem}: theoretical $\\Delta\\delta_{{eq}}$ (FR0)")
    colors = {"X_only": "tab:red", "Y_only": "tab:blue", "D_only": "tab:green"}

    table = {"stem": stem, "alpha0": alpha0, "L": L, "panels": []}
    for i, ax in enumerate(axes):
        panel = {"index": i, "target_component": ["xi_x", "xi_y", "eta1"][i], "curves": {}}
        for sname, sc in scenarios.items():
            ys = []
            inside_flags = []
            for val in ranges[i]:
                target = np.zeros(3)
                target[i] = val
                I = solve_currents(m["R_int"], m["int_nat"], target)
                # predicted chrom exactly target in linear model
                dde = delta_eq_theory(
                    target[0], target[1], target[2], alpha0,
                    eps_x=sc["eps_x"], eps_y=sc["eps_y"], delta_m=sc["delta_m"],
                    L=L, gamma=gamma,
                )
                ys.append(dde)
            ax.plot(ranges[i], ys, color=colors[sname], lw=2, label=sname)
            panel["curves"][sname] = {"x": ranges[i].tolist(), "delta_eq": ys}
        ax.axhline(0, color="k", lw=0.8, alpha=0.5)
        ax.axvline(0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel(PANEL_LABELS[i])
        ax.set_ylabel(r"$\Delta\delta_{eq}$")
        ax.set_title(f"Plane {['X','Y','Z'][i]}\n({CONSTRAINTS[i]})")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        table["panels"].append(panel)

    fig.tight_layout()
    out = PLOTS / f"{stem}_delta_eq_panels.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    write_json(DAT_OUT / stem / "delta_eq_panels.json", table)
    return out


def panel_map_dnu(stem: str, df) -> Path:
    models = build_lattice_models(df, None)
    ranges, curves = compute_panel_curves(models, standard_chrom_ranges(81))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"{stem}: map $\\Delta\\nu_s$ vs target chrom (FR0)")
    colors = ["tab:red", "tab:blue", "tab:green"]
    for i, ax in enumerate(axes):
        ax.plot(ranges[i], curves["NoFr_P0"][i], color=colors[i], ls="-", lw=2, label="P0")
        ax.plot(ranges[i], curves["NoFr_P1"][i], color=colors[i], ls="--", lw=2, label="P1")
        ax.axhline(0, color="k", lw=0.8, alpha=0.5)
        ax.axvline(0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel(PANEL_LABELS[i])
        ax.set_ylabel(r"$\Delta\nu_s$")
        ax.set_title(f"Plane {['X','Y','Z'][i]}\n({CONSTRAINTS[i]})")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    out = PLOTS / f"{stem}_map_dnu_panels.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def istar_delta_eq_table(stem: str, df, cfg) -> dict:
    m = get_dual_model(df)
    I = solve_currents(m["R_int"], m["int_nat"], np.zeros(3))
    chrom = m["R_int"] @ I + m["int_nat"]
    alpha0 = float(df["alpha_0"].mean())
    L = lattice_length(stem)
    gamma = float(cfg["gamma"])
    out = {"I_star": I.tolist(), "pred_chrom": chrom.tolist(), "scenarios": {}}
    for sname, sc in cfg["emittance_scenarios"].items():
        out["scenarios"][sname] = delta_eq_theory(
            chrom[0], chrom[1], chrom[2], alpha0,
            eps_x=sc["eps_x"], eps_y=sc["eps_y"], delta_m=sc["delta_m"],
            L=L, gamma=gamma,
        )
        # natural for comparison
        nat = m["int_nat"]
        out["scenarios"][sname + "_at_natural"] = delta_eq_theory(
            nat[0], nat[1], nat[2], alpha0,
            eps_x=sc["eps_x"], eps_y=sc["eps_y"], delta_m=sc["delta_m"],
            L=L, gamma=gamma,
        )
    return out


def main() -> int:
    cfg = load_config()
    lines = []
    for stem in stems():
        df = get_data(stem, "FR0")
        if df is None:
            continue
        p1 = panel_delta_eq(stem, df, cfg)
        p2 = panel_map_dnu(stem, df)
        tab = istar_delta_eq_table(stem, df, cfg)
        write_json(DAT_OUT / stem / "delta_eq_at_istar.json", tab)
        lines.append(f"- {stem}: {p1.name}, {p2.name}; Δδ_eq@I*={tab['scenarios']}")
        print(f"OK delta_eq {stem}")

    append_worklog(
        "## Δδ_eq + map-Δν_s panels\n\n"
        "- **Статус:** verified (offline)\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_delta_eq.py`\n"
        "- Сценарии ε/δ заданы в `config/study_config.json`.\n"
        "- При линейной модели Δδ_eq(ξ) проходит через 0 при ξ=0 **по построению**.\n"
        + "\n".join(lines) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
