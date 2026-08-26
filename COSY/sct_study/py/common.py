"""Shared paths and helpers for COSY/sct_study."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]
COSY = REPO / "COSY"
COSY_SRC = COSY / "src"
ANALYSIS = COSY / "analysis"
STRUCTURES = COSY / "structures"
DAT_SRC = COSY_SRC / "dat"
DAT_OUT = STUDY / "dat"
PLOTS = STUDY / "plots"
FOX = STUDY / "fox"
CONFIG_PATH = STUDY / "config" / "study_config.json"

sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(COSY_SRC / "run"))

from map_lat_lib import (  # noqa: E402
    get_data,
    get_dual_model,
    solve_currents,
    verify_chrom_model,
    verify_chrom_to_spin_link,
    audit_mapping_lattice,
    compute_panel_curves,
    build_lattice_models,
    standard_chrom_ranges,
)


def load_config() -> Dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def stems() -> List[str]:
    return list(load_config()["stems"])


def mapping_presentation_stems() -> List[str]:
    cfg = load_config()
    return list(cfg.get("mapping_presentation_stems", cfg["stems"]))


def lattice_fox_path(stem: str) -> Path:
    """Resolve lattice .fox from structures/ or COSY/src/."""
    candidates = [
        STRUCTURES / stem / f"{stem}.fox",
        COSY_SRC / f"{stem}.fox",
    ]
    for p in candidates:
        if p.is_file():
            return p
    raise FileNotFoundError(f"Lattice fox not found for {stem}")


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=_json_default), encoding="utf-8")


def _json_default(o: Any):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def append_worklog(section: str) -> None:
    log = STUDY / "WORKLOG.md"
    with log.open("a", encoding="utf-8") as f:
        f.write("\n" + section.rstrip() + "\n")


def update_manifest_job(job: Dict[str, Any]) -> None:
    man_path = STUDY / "manifest.json"
    man = json.loads(man_path.read_text(encoding="utf-8"))
    man.setdefault("jobs", []).append(job)
    man_path.write_text(json.dumps(man, indent=2), encoding="utf-8")


def lattice_length(stem: str) -> float:
    return float(load_config()["lattice_lengths_m"][stem])


def gamma_beta(cfg: Optional[Dict[str, Any]] = None):
    cfg = cfg or load_config()
    g = float(cfg["gamma"])
    b = float(np.sqrt(1.0 - 1.0 / g**2))
    return g, b


def delta_eq_theory(
    xi_x: float,
    xi_y: float,
    eta1: float,
    alpha0: float,
    *,
    eps_x: float,
    eps_y: float,
    delta_m: float,
    L: float,
    gamma: float,
) -> float:
    """Theoretical equilibrium relative-momentum shift Δδ_eq (Senichev et al.).

    This is *not* the COSY TRPRAY coordinate D and *not* the tracking proxy
    mean_D_offset = ⟨D⟩_i − ⟨D⟩_ref.
    """
    denom = gamma**2 * alpha0 - 1.0
    if abs(denom) < 1e-14:
        return float("nan")
    pre = gamma**2 / denom
    bracket = 0.5 * delta_m**2 * eta1 - (np.pi / L) * (eps_x * xi_x + eps_y * xi_y)
    return float(pre * bracket)


def mapping_grid_step(stem: str, mode: str = "FR0") -> np.ndarray:
    """Return FR0/FR3 current step vector (SGx1, SGx2, SGy1) for a stem."""
    from run_mapping import get_mapping_steps  # local import; ANALYSIS on path

    step_x, step_y = get_mapping_steps(stem, mode)
    # Scan uses SGx1,SGx2 along step_x and SGy1 along step_y (signed).
    return np.array([abs(step_x), abs(step_x), abs(step_y)], dtype=float)


def currents_in_box(
    I: np.ndarray,
    df,
    *,
    stem: Optional[str] = None,
    mode: str = "FR0",
) -> Dict[str, Any]:
    pts = df[["SGx1", "SGx2", "SGy1"]].values
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    inside = bool(np.all((I >= mins - 1e-15) & (I <= maxs + 1e-15)))
    dmin = float(np.linalg.norm(pts - I, axis=1).min())
    out: Dict[str, Any] = {
        "inside_axis_aligned_box": inside,
        "box_min": mins.tolist(),
        "box_max": maxs.tolist(),
        "min_distance_to_grid_point": dmin,
        "note": (
            "FR0 box = measured mapping domain, not a physical limit. "
            "Compare stems via min_distance_normalized, not raw d_min."
        ),
    }
    if stem is not None:
        step = mapping_grid_step(stem, mode)
        step_norm = float(np.linalg.norm(step))
        out["grid_step"] = step.tolist()
        out["grid_step_norm"] = step_norm
        out["min_distance_normalized"] = float(dmin / step_norm) if step_norm > 0 else float("nan")
    return out


def model_summary(df, stem: Optional[str] = None, mode: str = "FR0") -> Dict[str, Any]:
    m = get_dual_model(df, suffix="")
    m1 = get_dual_model(df, suffix="_1")
    Istar = solve_currents(m["R_int"], m["int_nat"], np.zeros(3))
    svals = np.linalg.svd(m["R_int"], compute_uv=False)
    cond = float(svals.max() / max(svals.min(), 1e-30))
    pred_spin = m["R_spin"] @ Istar + m["spin_nat"]
    pred_chrom = m["R_int"] @ Istar + m["int_nat"]
    box = currents_in_box(Istar, df, stem=stem, mode=mode)
    R = m["R_int"]
    # column correlations of response
    cols = R / np.linalg.norm(R, axis=0, keepdims=True)
    corr = cols.T @ cols
    return {
        "R_int": m["R_int"].tolist(),
        "int_nat": m["int_nat"].tolist(),
        "R_spin": m["R_spin"].tolist(),
        "spin_nat": m["spin_nat"].tolist(),
        "R_spin_P1": m1["R_spin"].tolist(),
        "spin_nat_P1": m1["spin_nat"].tolist(),
        "singular_values": svals.tolist(),
        "condition_number": cond,
        "column_correlation": corr.tolist(),
        "I_star": Istar.tolist(),
        "I_star_labels": ["SGx1", "SGx2", "SGy1"],
        "pred_chrom_at_Istar": pred_chrom.tolist(),
        "pred_dnu_s_P0_at_Istar": pred_spin.tolist(),
        "pred_dnu_s_P1_at_Istar": (m1["R_spin"] @ Istar + m1["spin_nat"]).tolist(),
        "grid_relation": box,
        "chrom_verify": verify_chrom_model(df),
        "spin_link_verify_P0": verify_chrom_to_spin_link(df, suffix=""),
        "spin_link_verify_P1": verify_chrom_to_spin_link(df, suffix="_1"),
        "alpha0_mean": float(df["alpha_0"].mean()),
        "alpha1_mean": float(df["alpha_1"].mean()),
        "eta1_mean": float(df["eta_1"].mean()),
        "n_points": int(len(df)),
        "SGy2_unique": int(df["SGy2"].nunique()),
    }
