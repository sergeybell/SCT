#!/usr/bin/env python3
"""Compute I* and offline predictions; prepare zero_point.json for COSY validation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DAT_OUT,
    append_worklog,
    get_data,
    get_dual_model,
    load_config,
    model_summary,
    solve_currents,
    stems,
    write_json,
)


def control_points(df, model) -> dict:
    """Pick in-grid currents with large |xi| for tracking comparisons."""
    # Use actual grid corners / mid extremes
    pts = []
    # natural
    pts.append({"tag": "natural", "I": [0.0, 0.0, 0.0], "note": "sextupoles off"})
    Istar = solve_currents(model["R_int"], model["int_nat"], np.zeros(3))
    pts.append({"tag": "Istar", "I": Istar.tolist(), "note": "target xi=eta1=0 (may extrapolate)"})

    # choose measured points maximizing |xi_x|, |xi_y|, |eta1| separately while keeping others moderate
    for col, tag in [("chrom_x", "ctrl_xi_x"), ("chrom_y", "ctrl_xi_y"), ("eta_1", "ctrl_eta1")]:
        idx = df[col].abs().idxmax()
        row = df.loc[idx]
        pts.append(
            {
                "tag": tag,
                "I": [float(row.SGx1), float(row.SGx2), float(row.SGy1)],
                "measured_chrom": [float(row.chrom_x), float(row.chrom_y), float(row.eta_1)],
                "note": f"grid point with max |{col}|",
            }
        )
    # annotate predicted chrom/spin for each
    for p in pts:
        I = np.array(p["I"])
        p["pred_chrom"] = (model["R_int"] @ I + model["int_nat"]).tolist()
        p["pred_dnu_s_P0"] = (model["R_spin"] @ I + model["spin_nat"]).tolist()
    return {"points": pts}


def main() -> int:
    rows = []
    for stem in stems():
        df = get_data(stem, "FR0")
        if df is None:
            continue
        summary = model_summary(df)
        m = get_dual_model(df)
        ctrls = control_points(df, m)
        out = {
            "stem": stem,
            "mode": "FR0",
            "I_star": summary["I_star"],
            "pred_chrom_at_Istar": summary["pred_chrom_at_Istar"],
            "pred_dnu_s_P0_at_Istar": summary["pred_dnu_s_P0_at_Istar"],
            "pred_dnu_s_P1_at_Istar": summary["pred_dnu_s_P1_at_Istar"],
            "condition_number": summary["condition_number"],
            "grid_relation": summary["grid_relation"],
            "alpha0_mean": summary["alpha0_mean"],
            "control_points": ctrls,
            "cosy_validation": {
                "status": "prepared",
                "fox": f"fox/validate_istar_{stem}.fox",
                "expected_dat": f"dat/{stem}/validate_istar.dat",
            },
        }
        write_json(DAT_OUT / stem / "zero_point.json", out)
        rows.append(
            f"| {stem} | {summary['I_star']} | {summary['pred_chrom_at_Istar']} | "
            f"{summary['pred_dnu_s_P0_at_Istar']} | {summary['grid_relation']['inside_axis_aligned_box']} |"
        )
        print(f"OK zero_point {stem}")

    append_worklog(
        "## Zero point I*\n\n"
        "- **Статус:** prepared (offline prediction); COSY direct check via validate_istar FOX\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_zero_point.py`\n"
        "| stem | I* | pred chrom | pred Δνs P0 | inside box |\n|------|----|------------|-------------|------------|\n"
        + "\n".join(rows)
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
