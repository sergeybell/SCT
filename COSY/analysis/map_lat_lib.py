"""Mapping lattice analysis utilities (shared by Map_magnetic.ipynb and plot_magnetic_map.py)."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = REPO / "COSY" / "src" / "dat"

COLS_INT = [
    "SGx1", "SGx2", "SGy1", "SGy2", "chrom_x", "chrom_y", "MU_ref",
    "quad_Kx", "quad_Ky", "quad_Kz", "alpha_0", "alpha_1", "eta_1",
]
# Delta spin tune from mapping.fox: MU_N_ARR(probe) - CONS(MU)
COLS_S0 = ["SGx1", "SGx2", "SGy1", "SGy2", "dnu_s_ref", "dnu_s_x", "dnu_s_y", "dnu_s_z"]
COLS_S1 = ["SGx1", "SGx2", "SGy1", "SGy2", "dnu_s_ref_1", "dnu_s_x_1", "dnu_s_y_1", "dnu_s_z_1"]
COLS_S0_FILE = ["SGx1", "SGx2", "SGy1", "SGy2", "spin_ref", "spin_x", "spin_y", "spin_z"]
COLS_S1_FILE = ["SGx1", "SGx2", "SGy1", "SGy2", "spin_ref_1", "spin_x_1", "spin_y_1", "spin_z_1"]

MAGNETIC_LATTICES = [f"magnetic_{n}" for n in range(2, 9)]

PANEL_CHROM_LABELS = [r"$\xi_x$", r"$\xi_y$", r"$\eta_1$"]
PANEL_CONSTRAINTS = [
    r"$\xi_y=0,\ \eta_1=0$",
    r"$\xi_x=0,\ \eta_1=0$",
    r"$\xi_x=0,\ \xi_y=0$",
]
PLANE_NAMES = ["x", "y", "z"]
PANEL_COLORS = ["red", "blue", "green"]
STEPS_DEFAULT = 100


def standard_chrom_ranges(steps: int = STEPS_DEFAULT) -> List[np.ndarray]:
    """Fixed target ranges as in Map_lat.ipynb: xi in [-10,10], eta_1 in [-1,1]."""
    return [
        np.linspace(-10, 10, steps),
        np.linspace(-10, 10, steps),
        np.linspace(-1, 1, steps),
    ]


def get_data(lattice: str, mode: str, root_path: Optional[Path] = None) -> Optional[pd.DataFrame]:
    base = Path(root_path or DEFAULT_ROOT) / lattice
    files = {
        "int": base / f"integrals_{mode}.dat",
        "s0": base / f"particle_spin_tune_{mode}.dat",
        "s1": base / f"particle_spin_tune_1_{mode}.dat",
    }
    if not all(p.is_file() for p in files.values()):
        print(f"Warning: missing files in {base} for mode {mode}")
        return None

    dfs = {}
    for key, path in files.items():
        file_cols = COLS_INT if key == "int" else (COLS_S0_FILE if key == "s0" else COLS_S1_FILE)
        canon = COLS_INT if key == "int" else (COLS_S0 if key == "s0" else COLS_S1)
        tmp = pd.read_csv(path, sep=r"\s+", names=file_cols, engine="python")
        tmp.columns = canon
        for col in tmp.columns:
            tmp[col] = pd.to_numeric(tmp[col].astype(str).str.replace("D", "E"), errors="coerce")
        dfs[key] = tmp.dropna()

    m = pd.merge(dfs["int"], dfs["s0"], on=["SGx1", "SGx2", "SGy1", "SGy2"])
    return pd.merge(m, dfs["s1"], on=["SGx1", "SGx2", "SGy1", "SGy2"])


def get_dual_model(df: pd.DataFrame, suffix: str = "") -> Dict[str, np.ndarray]:
    """Linear model: currents [SGx1, SGx2, SGy1] -> chrom [xi_x, xi_y, eta_1] and Delta nu_s probes."""
    X_design = np.column_stack([df["SGx1"], df["SGx2"], df["SGy1"], np.ones(len(df))])
    if suffix == "":
        y_dnu = df[["dnu_s_x", "dnu_s_y", "dnu_s_z"]].values
    else:
        y_dnu = df[[f"dnu_s_x{suffix}", f"dnu_s_y{suffix}", f"dnu_s_z{suffix}"]].values
    beta_s, _, _, _ = np.linalg.lstsq(X_design, y_dnu, rcond=None)
    y_chrom = df[["chrom_x", "chrom_y", "eta_1"]].values
    beta_c, _, _, _ = np.linalg.lstsq(X_design, y_chrom, rcond=None)
    return {
        "R_dnu_s": beta_s[:3, :].T,
        "dnu_s_nat": beta_s[3, :],
        "R_int": beta_c[:3, :].T,
        "int_nat": beta_c[3, :],
        "R_spin": beta_s[:3, :].T,
        "spin_nat": beta_s[3, :],
    }


def solve_currents(R_int: np.ndarray, int_nat: np.ndarray, target_chrom: np.ndarray) -> np.ndarray:
    rhs = target_chrom - int_nat
    try:
        return np.linalg.solve(R_int, rhs)
    except np.linalg.LinAlgError:
        I, _, _, _ = np.linalg.lstsq(R_int, rhs, rcond=None)
        return I


def build_lattice_models(
    df_no_fr: pd.DataFrame,
    df_with_fr: Optional[pd.DataFrame],
) -> Dict[str, Dict[str, np.ndarray]]:
    models = {
        "NoFr_P0": get_dual_model(df_no_fr, suffix=""),
        "NoFr_P1": get_dual_model(df_no_fr, suffix="_1"),
    }
    if df_with_fr is not None:
        models["WithFr_P0"] = get_dual_model(df_with_fr, suffix="")
        models["WithFr_P1"] = get_dual_model(df_with_fr, suffix="_1")
    return models


def compute_panel_curves(
    lattice_models: Dict[str, Dict[str, np.ndarray]],
    ranges: Optional[List[np.ndarray]] = None,
) -> Tuple[List[np.ndarray], Dict[str, List[List[float]]]]:
    """For each panel i, compute Delta nu_s curves vs target chrom component i."""
    ranges = ranges or standard_chrom_ranges()
    all_res: Dict[str, List[List[float]]] = {k: [] for k in lattice_models}

    for i in range(3):
        panel_res = {state: [] for state in lattice_models}
        m_no = lattice_models["NoFr_P0"]
        I_no = None
        I_with = None
        if "WithFr_P0" in lattice_models:
            m_with = lattice_models["WithFr_P0"]

        for val in ranges[i]:
            target = np.zeros(3)
            target[i] = val
            I_no = solve_currents(m_no["R_int"], m_no["int_nat"], target)
            if "WithFr_P0" in lattice_models:
                I_with = solve_currents(m_with["R_int"], m_with["int_nat"], target)
            for state, model in lattice_models.items():
                I = I_no if state.startswith("NoFr") else I_with
                dnu = float(np.dot(model["R_dnu_s"][i], I) + model["dnu_s_nat"][i])
                panel_res[state].append(dnu)

        for state in lattice_models:
            all_res[state].append(panel_res[state])

    return ranges, all_res
