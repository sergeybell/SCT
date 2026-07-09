"""Mapping lattice analysis utilities (shared by Map_magnetic.ipynb and plot_magnetic_map.py)."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = REPO / "COSY" / "src" / "dat"
STRUCTURES_ROOT = REPO / "COSY" / "structures"

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
CHROM_ABS_LIMIT = 20.0

# FR0: 5^3 grid (indices 0..4); FR3: 3^3 grid (indices 0..2)
GRID_MAX_INDEX = {"FR0": 4, "FR3": 2}


def scan_current_columns(df: pd.DataFrame) -> List[str]:
    """Independent sextupole gradients in the linear model (3 families)."""
    return ["SGx1", "SGx2", "SGy1"]


def _design_matrix(df: pd.DataFrame) -> np.ndarray:
    cols = scan_current_columns(df)
    parts = [df[c].values for c in cols] + [np.ones(len(df))]
    return np.column_stack(parts)


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
    """Linear model: sextupole currents -> chrom [xi_x, xi_y, eta_1] and Delta nu_s."""
    cols = scan_current_columns(df)
    X_design = _design_matrix(df)
    if suffix == "":
        y_dnu = df[["dnu_s_x", "dnu_s_y", "dnu_s_z"]].values
    else:
        y_dnu = df[[f"dnu_s_x{suffix}", f"dnu_s_y{suffix}", f"dnu_s_z{suffix}"]].values
    beta_s, _, _, _ = np.linalg.lstsq(X_design, y_dnu, rcond=None)
    y_chrom = df[["chrom_x", "chrom_y", "eta_1"]].values
    beta_c, _, _, _ = np.linalg.lstsq(X_design, y_chrom, rcond=None)
    n = len(cols)
    return {
        "scan_cols": cols,
        "R_dnu_s": beta_s[:n, :].T,
        "dnu_s_nat": beta_s[n, :],
        "R_int": beta_c[:n, :].T,
        "int_nat": beta_c[n, :],
        "R_spin": beta_s[:n, :].T,
        "spin_nat": beta_s[n, :],
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


def _r2_score(actual: np.ndarray, predicted: np.ndarray) -> float:
    ss_res = np.sum((actual - predicted) ** 2)
    ss_tot = np.sum((actual - np.mean(actual)) ** 2)
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0
    return float(1.0 - ss_res / ss_tot)


def predict_chrom(df: pd.DataFrame, model: Optional[Dict[str, np.ndarray]] = None) -> np.ndarray:
    """Predict chromaticity from sextupole currents using the linear model."""
    m = model or get_dual_model(df, suffix="")
    cols = m.get("scan_cols", scan_current_columns(df))
    currents = df[cols].values
    return np.array([np.dot(m["R_int"], i) + m["int_nat"] for i in currents])


def predict_spin_via_chrom(
    df: pd.DataFrame,
    suffix: str = "",
    model: Optional[Dict[str, np.ndarray]] = None,
) -> np.ndarray:
    """Predict Delta nu_s by inverting chrom model: xi -> I -> spin (Map_lat cell 13)."""
    m = model or get_dual_model(df, suffix=suffix)
    if suffix == "":
        spin_cols = ["dnu_s_x", "dnu_s_y", "dnu_s_z"]
    else:
        spin_cols = [f"dnu_s_x{suffix}", f"dnu_s_y{suffix}", f"dnu_s_z{suffix}"]
    xi_actual = df[["chrom_x", "chrom_y", "eta_1"]].values
    predicted = []
    for xi in xi_actual:
        i_calc = solve_currents(m["R_int"], m["int_nat"], xi)
        predicted.append(np.dot(m["R_spin"], i_calc) + m["spin_nat"])
    return np.array(predicted), df[spin_cols].values


def chrom_range_report(
    df: pd.DataFrame,
    *,
    limit: float = CHROM_ABS_LIMIT,
) -> Dict[str, Union[float, bool]]:
    """Summarize chromaticity ranges and flag transverse excursions beyond limit."""
    report = {
        "chrom_x_min": float(df["chrom_x"].min()),
        "chrom_x_max": float(df["chrom_x"].max()),
        "chrom_y_min": float(df["chrom_y"].min()),
        "chrom_y_max": float(df["chrom_y"].max()),
        "eta_1_min": float(df["eta_1"].min()),
        "eta_1_max": float(df["eta_1"].max()),
        "n_points": len(df),
    }
    report["xi_in_range"] = (
        report["chrom_x_min"] >= -limit
        and report["chrom_x_max"] <= limit
        and report["chrom_y_min"] >= -limit
        and report["chrom_y_max"] <= limit
    )
    return report


def verify_chrom_model(df: pd.DataFrame) -> Dict[str, float]:
    """Check linear chrom model: measured vs predicted from currents."""
    actual = df[["chrom_x", "chrom_y", "eta_1"]].values
    predicted = predict_chrom(df)
    names = ["chrom_x", "chrom_y", "eta_1"]
    out: Dict[str, float] = {}
    for j, name in enumerate(names):
        out[f"r2_{name}"] = _r2_score(actual[:, j], predicted[:, j])
        out[f"max_res_{name}"] = float(np.max(np.abs(actual[:, j] - predicted[:, j])))
    return out


def verify_chrom_to_spin_link(
    df: pd.DataFrame,
    suffix: str = "",
) -> Dict[str, float]:
    """Check xi -> spin chain (Map_lat.ipynb cell 13)."""
    predicted, actual = predict_spin_via_chrom(df, suffix=suffix)
    planes = ["x", "y", "z"]
    out: Dict[str, float] = {}
    for j, plane in enumerate(planes):
        out[f"r2_{plane}"] = _r2_score(actual[:, j], predicted[:, j])
        out[f"max_res_{plane}"] = float(np.max(np.abs(actual[:, j] - predicted[:, j])))
    return out


def estimate_mapping_steps(
    df: pd.DataFrame,
    mode: str,
    *,
    limit: float = CHROM_ABS_LIMIT,
    step_ratio: float = -1.5,
) -> Tuple[float, float]:
    """Suggest step_x/step_y so grid corner stays within |xi| <= limit."""
    m = get_dual_model(df, suffix="")
    max_idx = GRID_MAX_INDEX.get(mode, 4)
    sx = float(df["SGx2"].replace(0, np.nan).dropna().diff().abs().min())
    if not np.isfinite(sx) or sx == 0:
        sx = 0.02
    corner = np.array([max_idx * sx, max_idx * sx, max_idx * step_ratio * sx])
    chrom_corner = np.dot(m["R_int"], corner) + m["int_nat"]
    scale = max(np.max(np.abs(chrom_corner[:2])), 1e-9) / limit
    step_x = sx / scale
    step_y = step_ratio * step_x
    return round(step_x, 6), round(step_y, 6)


def plot_verify_chrom_to_spin(
    lattice: str,
    mode: str,
    df: pd.DataFrame,
    *,
    suffix: str = "",
    show: bool = False,
    save_dir: Optional[Path] = None,
):
    """Scatter measured vs predicted Delta nu_s (chrom-mediated link)."""
    import matplotlib.pyplot as plt

    predicted, actual = predict_spin_via_chrom(df, suffix=suffix)
    particle = "P0" if suffix == "" else "P1"
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(
        f"Link verification: $\\xi \\to \\Delta\\nu_s$ — {lattice} {mode} ({particle})",
        fontsize=14,
    )
    plane_labels = ["X", "Y", "Z"]
    for i, ax in enumerate(axes):
        true_s = actual[:, i]
        pred_s = predicted[:, i]
        r2 = _r2_score(true_s, pred_s)
        ax.scatter(true_s, pred_s, color="purple", alpha=0.5, label="Actual link")
        lims = [min(true_s.min(), pred_s.min()), max(true_s.max(), pred_s.max())]
        ax.plot(lims, lims, "r--", label="Perfect model")
        ax.set_title(f"Plane {plane_labels[i]}\n$R^2 = {r2:.8f}$")
        ax.set_xlabel("Measured $\\Delta\\nu_s$ (COSY)")
        ax.set_ylabel("Predicted $\\Delta\\nu_s$ (via $\\xi$)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        out = save_dir / f"{lattice}_verify_spin_link_{mode}_{particle}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"OK: {out}")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def audit_mapping_lattice(lattice: str, root_path: Optional[Path] = None) -> Dict[str, Union[int, bool, str, float, Tuple[float, float]]]:
    """Summarize 3-family sextupole wiring vs mapping scan for one magnetic structure."""
    fox_path = STRUCTURES_ROOT / lattice / f"{lattice}.fox"
    if not fox_path.is_file():
        return {"error": f"missing {fox_path}"}

    text = fox_path.read_text(encoding="utf-8")
    sf1 = text.count("MH L_SF1")
    sf2 = text.count("MH L_SF2")
    sd = text.count("MH L_SD SD")
    wired = {
        "SF1_from_SEXTGx1": "SF1 := SEXTGx1" in text,
        "SF2_from_SEXTGx2": "SF2 := SEXTGx2" in text,
        "SD_from_SEXTGy1": "SD := SEXTGy1" in text,
        "SD1_present": "SD1 := SEXTGy2" in text or "MH L_SD1" in text,
    }
    if wired["SD1_present"]:
        print(f"WARNING {lattice}: unexpected 4th sextupole family SD1 — should be 3 families only")

    out: Dict[str, Union[int, bool, str, float, Tuple[float, float]]] = {
        "lattice": lattice,
        "physical_SF1": sf1,
        "physical_SF2": sf2,
        "physical_SD": sd,
        "physical_sextupoles_total": sf1 + sf2 + sd,
        "independent_powered_families": 3,
        "scan_parameters": "SGx1, SGx2, SGy1",
        "FR0_grid_points": 125,
        "FR3_grid_points": 27,
        **wired,
    }

    root = Path(root_path or DEFAULT_ROOT)
    for mode in ("FR0", "FR3"):
        df = get_data(lattice, mode, root)
        if df is None:
            out[f"{mode}_status"] = "missing dat"
            continue
        rep = chrom_range_report(df)
        out[f"{mode}_n_points"] = rep["n_points"]
        out[f"{mode}_SGy2_varies"] = df["SGy2"].nunique() > 1
        out[f"{mode}_xi_in_range"] = rep["xi_in_range"]
        out[f"{mode}_xi_x_max"] = max(abs(rep["chrom_x_min"]), abs(rep["chrom_x_max"]))
        out[f"{mode}_xi_y_max"] = max(abs(rep["chrom_y_min"]), abs(rep["chrom_y_max"]))
        est = estimate_mapping_steps(df, mode)
        out[f"{mode}_suggested_steps"] = est

    return out


def print_mapping_audit(lattice: str, root_path: Optional[Path] = None) -> None:
    """Human-readable audit of mapping setup vs data."""
    a = audit_mapping_lattice(lattice, root_path)
    if "error" in a:
        print(a["error"])
        return
    print(f"\n=== Mapping audit: {lattice} ===")
    print(
        f"Physical sextupoles: {a['physical_sextupoles_total']} "
        f"(SF1={a['physical_SF1']}, SF2={a['physical_SF2']}, SD={a['physical_SD']})"
    )
    print(f"Powered families: {a['independent_powered_families']}  |  Scan varies: {a['scan_parameters']}")
    print(f"Wiring: SF1<-SEXTGx1={a['SF1_from_SEXTGx1']}, SF2<-SEXTGx2={a['SF2_from_SEXTGx2']}, "
          f"SD<-SEXTGy1={a['SD_from_SEXTGy1']}, SD1_present={a['SD1_present']}")
    print(f"Grid: FR0={a['FR0_grid_points']} pts, FR3={a['FR3_grid_points']} pts")
    for mode in ("FR0", "FR3"):
        key = f"{mode}_status"
        if key in a:
            print(f"  {mode}: {a[key]}")
            continue
        print(
            f"  {mode}: n={a[f'{mode}_n_points']}, SGy2 varies={a[f'{mode}_SGy2_varies']}, "
            f"|xi| max ~{max(a[f'{mode}_xi_x_max'], a[f'{mode}_xi_y_max']):.1f}, "
            f"in_range={a[f'{mode}_xi_in_range']}, suggested steps={a[f'{mode}_suggested_steps']}"
        )
