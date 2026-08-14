#!/usr/bin/env python3
"""OptiM NuX/NuY polar vs Twiss ψ and geometric 2π Q_full s/C (magnetic_2).

Nu is betatron oscillations (turns), not an angle: φ = 2π Nu.
Amplitudes: |I*| with lattice signs SF1>0, SF2>0, SD<0 (I* itself may flip SF2).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_phase_advance import (  # noqa: E402
    default_tunes,
    resultant_vector,
    signed_polar,
)
from common import (  # noqa: E402
    DAT_OUT,
    PLOTS,
    REPO,
    append_worklog,
    lattice_length,
    write_json,
)

STEM = "magnetic_2"
OPTIM_TABLE = REPO / "OptiM" / "magnetic" / "sext_phase_magnetic_2.txt"
COLORS = {"SF1": "tab:red", "SF2": "tab:orange", "SD": "tab:blue"}


def wrap_turns(d: float) -> float:
    return float((d + 0.5) % 1.0 - 0.5)


def family_amplitudes_from_istar(Istar: np.ndarray) -> Dict[str, float]:
    """Focusing SF1, SF2 > 0; defocusing SD < 0. Magnitude from |I*|."""
    return {
        "SF1": abs(float(Istar[0])),
        "SF2": abs(float(Istar[1])),
        "SD": -abs(float(Istar[2])),
    }


def load_istar(stem: str) -> np.ndarray:
    zp = DAT_OUT / stem / "zero_point.json"
    data = json.loads(zp.read_text(encoding="utf-8"))
    return np.asarray(data["I_star"], dtype=float)


def parse_optim_table(path: Path) -> Tuple[List[dict], Dict[str, float]]:
    """Parse OptiM Twiss dump. S[cm] and Nu are at element exit."""
    rows: List[dict] = []
    end: Optional[dict] = None
    text = path.read_text(encoding="utf-8")
    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split()
        if len(parts) < 13:
            continue
        rec = {
            "n": int(float(parts[0])),
            "name": parts[1],
            "s_cm": float(parts[2]),
            "beta_x_cm": float(parts[3]),
            "beta_y_cm": float(parts[5]),
            "nu_x": float(parts[11]),
            "nu_y": float(parts[12]),
        }
        rec["s_m"] = rec["s_cm"] / 100.0
        if rec["name"] == "END":
            end = rec
        rows.append(rec)
    if end is None:
        raise ValueError(f"no END row in {path}")
    return rows, {"Qx": end["nu_x"], "Qy": end["nu_y"], "C_m": end["s_m"]}


def extract_sextupoles(rows: List[dict], strengths: Dict[str, float]) -> List[dict]:
    out = []
    for rec in rows:
        fam = rec["name"]
        if fam not in strengths:
            continue
        e = dict(rec)
        e["family"] = fam
        e["K2"] = float(strengths[fam])
        e["s_center"] = e["s_m"]  # OptiM table is exit; label plots with that s
        e["mu_x_rad"] = 2.0 * np.pi * e["nu_x"]
        e["mu_y_rad"] = 2.0 * np.pi * e["nu_y"]
        out.append(e)
    return out


def load_twiss(stem: str) -> Optional[dict]:
    path = DAT_OUT / stem / "sext_phase_twiss.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def match_twiss(optim: List[dict], twiss_sexts: List[dict], tol: float = 2e-3) -> List[dict]:
    if len(optim) != len(twiss_sexts):
        raise ValueError(f"count mismatch: OptiM {len(optim)} vs Twiss {len(twiss_sexts)}")
    matched = []
    for o, t in zip(optim, twiss_sexts):
        if o["family"] != t["family"]:
            raise ValueError(f"family mismatch: {o['family']} vs {t['family']}")
        ds = abs(o["s_m"] - float(t["s_end"]))
        if ds > tol:
            raise ValueError(
                f"{o['family']} s mismatch: OptiM {o['s_m']:.6f} vs COSY s_end {t['s_end']:.6f} (Δ={ds:.4g})"
            )
        matched.append(t)
    return matched


def with_phase(sexts: List[dict], mux: List[float], muy: List[float]) -> List[dict]:
    out = []
    for el, ax, ay in zip(sexts, mux, muy):
        e = dict(el)
        e["mu_x_rad"] = float(ax)
        e["mu_y_rad"] = float(ay)
        out.append(e)
    return out


def plot_signed_polar(
    sexts: List[dict],
    plane: str,
    resultant: dict,
    title: str,
    outfile: Path,
) -> Path:
    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="polar")
    radii = []
    for el in sexts:
        fam = el["family"]
        mu = el["mu_x_rad"] if plane == "x" else el["mu_y_rad"]
        ang, r = signed_polar(el["K2"], mu)
        radii.append(r)
        ax.annotate(
            "",
            xy=(ang, r),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=COLORS.get(fam, "k"), lw=1.8),
        )
        ax.text(
            ang,
            r * 1.08,
            f"{fam}\ns={el['s_m']:.2f} m",
            fontsize=7,
            ha="center",
            va="bottom",
        )
    r_sum = float(resultant["magnitude"])
    a_sum = float(resultant["angle_rad"])
    radii.append(r_sum)
    ax.annotate(
        "",
        xy=(a_sum, r_sum),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color="black", lw=2.6, linestyle="--"),
    )
    ax.text(a_sum, max(r_sum, 1e-12) * 1.12, "Σ", fontsize=11, fontweight="bold", ha="center")
    rmax = max(radii) if radii else 1.0
    ax.set_rlim(0, rmax * 1.25 if rmax > 0 else 1.0)
    ax.set_title(title, pad=18)
    for fam, col in COLORS.items():
        ax.plot([], [], color=col, lw=2, label=fam)
    ax.plot([], [], color="black", lw=2.6, ls="--", label="resultant Σ")
    ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08), fontsize=8)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outfile


def plot_compare_row(
    panels: List[Tuple[str, List[dict], dict]],
    plane: str,
    outfile: Path,
) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6), subplot_kw={"projection": "polar"})
    rmax = 0.0
    drawn = []
    for ax, (title, sexts, res) in zip(axes, panels):
        radii = []
        for el in sexts:
            fam = el["family"]
            mu = el["mu_x_rad"] if plane == "x" else el["mu_y_rad"]
            ang, r = signed_polar(el["K2"], mu)
            radii.append(r)
            ax.annotate(
                "",
                xy=(ang, r),
                xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color=COLORS.get(fam, "k"), lw=1.6),
            )
        r_sum = float(res["magnitude"])
        a_sum = float(res["angle_rad"])
        radii.append(r_sum)
        ax.annotate(
            "",
            xy=(a_sum, r_sum),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2.2, linestyle="--"),
        )
        local = max(radii) if radii else 1.0
        rmax = max(rmax, local)
        drawn.append((ax, title, res))
    for ax, title, res in drawn:
        ax.set_rlim(0, rmax * 1.25 if rmax > 0 else 1.0)
        ax.set_title(
            f"{title}\n|Σ|={res['magnitude']:.4g}, ∠={res['angle_deg']:.1f}°",
            pad=14,
            fontsize=10,
        )
    for fam, col in COLORS.items():
        axes[0].plot([], [], color=col, lw=2, label=fam)
    axes[0].plot([], [], color="black", lw=2.2, ls="--", label="Σ")
    axes[0].legend(loc="upper right", bbox_to_anchor=(1.45, 1.12), fontsize=8)
    fig.suptitle(f"{STEM}: A·e^{{iφ}}  ({plane})", y=1.04, fontsize=12)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outfile


def rms(vals: List[float]) -> float:
    arr = np.asarray(vals, dtype=float)
    return float(np.sqrt(np.mean(arr**2))) if len(arr) else float("nan")


def analyze() -> dict:
    Istar = load_istar(STEM)
    Istar_raw = {
        "SF1": float(Istar[0]),
        "SF2": float(Istar[1]),
        "SD": float(Istar[2]),
    }
    strengths = family_amplitudes_from_istar(Istar)
    rows, ring = parse_optim_table(OPTIM_TABLE)
    C = lattice_length(STEM)
    Qx_full = float(ring["Qx"])
    Qy_full = float(ring["Qy"])
    Q_int = int(np.floor(min(Qx_full, Qy_full)))
    Qx_frac, Qy_frac = default_tunes(STEM)
    Qx_cosy_full = Q_int + Qx_frac
    Qy_cosy_full = Q_int + Qy_frac

    optim = extract_sextupoles(rows, strengths)
    if not optim:
        raise RuntimeError("no SF1/SF2/SD rows in OptiM table")

    twiss_doc = load_twiss(STEM)
    twiss_matched = match_twiss(optim, twiss_doc["sextupoles"]) if twiss_doc else []

    rows_out = []
    dx_twiss, dy_twiss, dx_geom, dy_geom = [], [], [], []
    for i, o in enumerate(optim):
        s = o["s_m"]
        frac = s / C if C else 0.0
        geom_x = Qx_full * frac
        geom_y = Qy_full * frac
        frac_x = Qx_frac * frac
        frac_y = Qy_frac * frac
        rec = {
            "family": o["family"],
            "s_m": s,
            "K2": o["K2"],
            "nu_x": o["nu_x"],
            "nu_y": o["nu_y"],
            "geom_full_x_turns": geom_x,
            "geom_full_y_turns": geom_y,
            "geom_frac_x_turns": frac_x,
            "geom_frac_y_turns": frac_y,
            "d_nu_x_geom_full": wrap_turns(geom_x - o["nu_x"]),
            "d_nu_y_geom_full": wrap_turns(geom_y - o["nu_y"]),
        }
        if twiss_matched:
            t = twiss_matched[i]
            rec["psi_x_turns"] = float(t["psi_x_turns"])
            rec["psi_y_turns"] = float(t["psi_y_turns"])
            rec["s_end_cosy"] = float(t["s_end"])
            rec["d_nu_x_twiss"] = wrap_turns(t["psi_x_turns"] - o["nu_x"])
            rec["d_nu_y_twiss"] = wrap_turns(t["psi_y_turns"] - o["nu_y"])
            dx_twiss.append(rec["d_nu_x_twiss"])
            dy_twiss.append(rec["d_nu_y_twiss"])
        dx_geom.append(rec["d_nu_x_geom_full"])
        dy_geom.append(rec["d_nu_y_geom_full"])
        rows_out.append(rec)

    geom_sexts = with_phase(
        optim,
        [2 * np.pi * r["geom_full_x_turns"] for r in rows_out],
        [2 * np.pi * r["geom_full_y_turns"] for r in rows_out],
    )
    twiss_sexts = []
    if twiss_matched:
        twiss_sexts = with_phase(
            optim,
            [2 * np.pi * r["psi_x_turns"] for r in rows_out],
            [2 * np.pi * r["psi_y_turns"] for r in rows_out],
        )

    methods = {
        "optim": {
            "sexts": optim,
            "note": "φ=2π Nu from OptiM table (turns, element exit)",
        },
        "geom_full": {
            "sexts": geom_sexts,
            "note": "φ=2π Q_full s/C; Q_full = Nu(END)",
        },
    }
    if twiss_sexts:
        methods["twiss"] = {
            "sexts": twiss_sexts,
            "note": "φ=ψ=∫ds/β from COSY Twiss (element center)",
        }

    resultants = {}
    plots = {}
    for key, spec in methods.items():
        rx = resultant_vector(spec["sexts"], "x")
        ry = resultant_vector(spec["sexts"], "y")
        rx["note"] = spec["note"]
        ry["note"] = spec["note"]
        resultants[key] = {"x": rx, "y": ry}

    plots["optim_x"] = str(
        plot_signed_polar(
            optim,
            "x",
            resultants["optim"]["x"],
            f"{STEM}: OptiM  A e^{{i 2π Nu_x}}\n"
            f"|Σ|={resultants['optim']['x']['magnitude']:.4g},  "
            f"∠Σ={resultants['optim']['x']['angle_deg']:.1f}°",
            PLOTS / f"{STEM}_sext_phase_optim_x.png",
        )
    )
    plots["optim_y"] = str(
        plot_signed_polar(
            optim,
            "y",
            resultants["optim"]["y"],
            f"{STEM}: OptiM  A e^{{i 2π Nu_y}}\n"
            f"|Σ|={resultants['optim']['y']['magnitude']:.4g},  "
            f"∠Σ={resultants['optim']['y']['angle_deg']:.1f}°",
            PLOTS / f"{STEM}_sext_phase_optim_y.png",
        )
    )
    plots["geom_full_x"] = str(
        plot_signed_polar(
            geom_sexts,
            "x",
            resultants["geom_full"]["x"],
            f"{STEM}: geom  A e^{{i 2π Q_x s/C}}  (Q_x={Qx_full:.5f})\n"
            f"|Σ|={resultants['geom_full']['x']['magnitude']:.4g},  "
            f"∠Σ={resultants['geom_full']['x']['angle_deg']:.1f}°",
            PLOTS / f"{STEM}_sext_phase_geom_full_x.png",
        )
    )
    plots["geom_full_y"] = str(
        plot_signed_polar(
            geom_sexts,
            "y",
            resultants["geom_full"]["y"],
            f"{STEM}: geom  A e^{{i 2π Q_y s/C}}  (Q_y={Qy_full:.5f})\n"
            f"|Σ|={resultants['geom_full']['y']['magnitude']:.4g},  "
            f"∠Σ={resultants['geom_full']['y']['angle_deg']:.1f}°",
            PLOTS / f"{STEM}_sext_phase_geom_full_y.png",
        )
    )
    if twiss_sexts:
        for plane in ("x", "y"):
            panels = [
                ("OptiM Nu", optim, resultants["optim"][plane]),
                ("Twiss ψ", twiss_sexts, resultants["twiss"][plane]),
                ("Q_full·s/C", geom_sexts, resultants["geom_full"][plane]),
            ]
            plots[f"compare_{plane}"] = str(
                plot_compare_row(panels, plane, PLOTS / f"{STEM}_sext_phase_compare_{plane}.png")
            )

    out = {
        "stem": STEM,
        "source_table": str(OPTIM_TABLE),
        "n_sextupole_instances": len(optim),
        "C_m": C,
        "C_optim_end_m": ring["C_m"],
        "I_star_raw": Istar_raw,
        "amplitudes": strengths,
        "amplitude_note": "A = |I*| with lattice signs SF1>0, SF2>0, SD<0; I* solve may flip SF2",
        "Qx_full_optim": Qx_full,
        "Qy_full_optim": Qy_full,
        "Q_int": Q_int,
        "Qx_frac_cosy": Qx_frac,
        "Qy_frac_cosy": Qy_frac,
        "Qx_full_cosy": Qx_cosy_full,
        "Qy_full_cosy": Qy_cosy_full,
        "phase_note": "Nu is turns; polar angle φ=2π Nu. Legacy sext_phase.json uses Q_frac only.",
        "rms_dnu_x_twiss": rms(dx_twiss) if dx_twiss else None,
        "rms_dnu_y_twiss": rms(dy_twiss) if dy_twiss else None,
        "rms_dnu_x_geom_full": rms(dx_geom),
        "rms_dnu_y_geom_full": rms(dy_geom),
        "sextupoles": rows_out,
        "resultant": resultants,
        "plots": plots,
    }
    write_json(DAT_OUT / STEM / "sext_phase_optim.json", out)
    return out


def main() -> int:
    if not OPTIM_TABLE.is_file():
        print(f"SKIP: missing {OPTIM_TABLE}")
        return 0
    out = analyze()
    print(
        f"OK optim-phase {STEM}: {out['n_sextupole_instances']} instances; "
        f"A={out['amplitudes']}; "
        f"rms|dNu|_geom x/y={out['rms_dnu_x_geom_full']:.4g}/{out['rms_dnu_y_geom_full']:.4g}"
    )
    if out["rms_dnu_x_twiss"] is not None:
        print(
            f"  rms|dNu|_twiss x/y={out['rms_dnu_x_twiss']:.4g}/{out['rms_dnu_y_twiss']:.4g}"
        )
    append_worklog(
        "## OptiM Nu vs Twiss vs Q_full·s/C\n\n"
        f"- **Статус:** offline from `{OPTIM_TABLE.name}` + I* + existing Twiss JSON\n"
        f"- **Команда:** `python COSY/sct_study/py/analyze_optim_sext_phase.py`\n"
        "- A: |I*| with SF1>0, SF2>0, SD<0. φ=2π Nu (Nu in turns).\n"
        f"- {STEM}: {out['n_sextupole_instances']} instances; "
        f"rms|ΔNu|_geom={out['rms_dnu_x_geom_full']:.4g}/{out['rms_dnu_y_geom_full']:.4g}; "
        f"rms|ΔNu|_twiss={out.get('rms_dnu_x_twiss')}/{out.get('rms_dnu_y_twiss')}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
