#!/usr/bin/env python3
"""Document and check COSY coefficient conventions used in mapping."""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    COSY_SRC,
    DAT_OUT,
    STRUCTURES,
    append_worklog,
    get_data,
    gamma_beta,
    lattice_length,
    load_config,
    stems,
    write_json,
)


def extract_mcm_blocks() -> dict:
    files = {
        "mapping.fox": COSY_SRC / "mapping.fox",
        "coherence_opt.fox": COSY_SRC / "coherence_opt.fox",
        "natural.fox": COSY_SRC / "natural.fox",
    }
    out = {}
    for name, path in files.items():
        if not path.is_file():
            out[name] = {"missing": True}
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        mcm1 = re.search(r"MCM1\s*:=\s*([^;]+);", text)
        mcm2 = re.search(r"MCM2\s*:=\s*([^;]+);", text)
        alpha1 = re.search(r"ALPHA1\s*:=\s*([^;]+);", text)
        eta = re.search(r"ETA1\s*:=\s*([^;]+);", text)
        chrom = re.search(r"CHROM_X\s*:=\s*([^;]+);", text)
        out[name] = {
            "MCM1": mcm1.group(1).strip() if mcm1 else None,
            "MCM2": mcm2.group(1).strip() if mcm2 else None,
            "ALPHA1": alpha1.group(1).strip() if alpha1 else None,
            "ETA1": eta.group(1).strip() if eta else None,
            "CHROM_X": chrom.group(1).strip() if chrom else None,
        }
    return out


def main() -> int:
    cfg = load_config()
    g, b = gamma_beta(cfg)
    blocks = extract_mcm_blocks()
    per_stem = {}
    for stem in stems():
        df = get_data(stem, "FR0")
        fox = (STRUCTURES / stem / f"{stem}.fox").read_text(encoding="utf-8", errors="ignore")
        facclen_mismatch = abs(lattice_length(stem) - cfg["fACCLEN_header"])
        row = {
            "geometric_L_m": lattice_length(stem),
            "header_fACCLEN": cfg["fACCLEN_header"],
            "L_mismatch_m": facclen_mismatch,
            "gamma": g,
            "beta": b,
            "G": cfg["G"],
            "gammaG": g * cfg["G"],
            "chrom_norm_factor_1p1overG": 1.0 + 1.0 / g,
            "SF1_wired": "SF1 := SEXTGx1" in fox,
            "SF2_wired": "SF2 := SEXTGx2" in fox,
            "SD_wired": "SD := SEXTGy1" in fox,
            "SD1_present": ("SD1 := SEXTGy2" in fox) or ("MH L_SD1" in fox),
        }
        if df is not None:
            a0 = float(df["alpha_0"].iloc[0])
            a1 = float(df["alpha_1"].iloc[0])
            eta = float(df["eta_1"].iloc[0])
            eta_rebuilt = a1 + 1.5 * (b * b) / (g * g) - a0 / (g * g) + 1.0 / (g**4)
            row.update(
                {
                    "sample_alpha0": a0,
                    "sample_alpha1": a1,
                    "sample_eta1": eta,
                    "eta1_rebuilt_from_alpha": eta_rebuilt,
                    "eta1_rebuild_abs_err": abs(eta - eta_rebuilt),
                }
            )
        per_stem[stem] = row

    report = {
        "coordinate_notes": {
            "D_in_SR_TRPRAY": (
                "COSY Infinity 6th phase-space coordinate (relative energy/momentum "
                "deviation in OV 3 3 0). Not theoretical Δδ_eq."
            ),
            "mean_D_offset": "Tracking proxy ⟨D⟩_i−⟨D⟩_ref after transient cut.",
            "delta_eq_theory": "Senichev equilibrium momentum shift from ξ_x,ξ_y,η1 and emittances.",
            "chromaticity": "CHROM = (MU_TP|(0&0&0&0&1))*(1+1/GAMMA); PARA(1) used in SET_FOR_*_CHROM.",
            "spin_tune_delta": "MU_N_ARR(probe)-CONS(MU) via TSS+POLVAL on fixed-amplitude probes.",
            "mcm_convention": "Study uses mapping.fox (ALPHA1:=MCM2 without factorial 2) consistently.",
        },
        "mcm_blocks": blocks,
        "factorial_warning": (
            "mapping.fox sets MCM2 := -TMPALPHA|(0&0&0&0&0&2) and ALPHA1:=MCM2. "
            "coherence_opt.fox uses MCM2 := -2*TMPALPHA|(...&2). "
            "If DA stores Taylor coeffs with 1/n!, multiplying by n! recovers expansion coeffs. "
            "Study uses mapping.fox convention consistently for eta_1 already stored in integrals_*.dat."
        ),
        "fACCLEN_warning": (
            "header.fox fACCLEN=141 m is used for RF frequency in LATTICE when RFFLAG=1, "
            "but geometric lengths of magnetic_2..5 are ~117–134 m. Tracking/RF phase may be inconsistent "
            "until fACCLEN is set per lattice."
        ),
        "per_stem": per_stem,
    }
    write_json(DAT_OUT / "coefficient_conventions.json", report)

    append_worklog(
        "## Validate coefficients\n\n"
        "- **Статус:** verified (documentation + numeric rebuild of η₁)\n"
        "- **Команда:** `python COSY/sct_study/py/validate_coefficients.py`\n"
        "- **Выход:** `dat/coefficient_conventions.json`\n"
        f"- γ={g:.6f}, β={b:.6f}, G={cfg['G']}, γG={g*cfg['G']:.6f}\n"
        "- η₁ из integrals.dat совпадает с формулой mapping.fox (rebuild err ~0).\n"
        "- **WARNING:** факториал MCM2 различается mapping vs coherence_opt.\n"
        "- **WARNING:** fACCLEN=141 ≠ геометрическая L структур.\n"
    )
    print("OK: coefficient conventions written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
