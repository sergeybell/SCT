#!/usr/bin/env python3
"""Twiss-based betatron phase at each sextupole instance vs geometric s/C·Q proxy."""
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
    parse_lengths,
    parse_sext_strengths_from_istar,
    plot_polar,
    resultant_vector,
    walk_lattice,
)
from common import (  # noqa: E402
    DAT_OUT,
    DAT_SRC,
    PLOTS,
    STRUCTURES,
    append_worklog,
    get_data,
    get_dual_model,
    solve_currents,
    stems,
    write_json,
)


def read_beta(path: Path) -> Dict[int, float]:
    out: Dict[int, float] = {}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.upper().startswith("BETA") or line.upper().startswith("DISP"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            beta = float(parts[0].replace("D", "E"))
            idx = int(float(parts[1].replace("D", "E")))
        except ValueError:
            continue
        out[idx] = beta
    return out


def element_index_map(elements: List[dict]) -> List[dict]:
    """Assign COSY BETS element indices: 0 = start, j after j-th element."""
    out = []
    for j, el in enumerate(elements, start=1):
        e = dict(el)
        e["ele_index"] = j
        out.append(e)
    return out


def integrate_phase(elements: List[dict], beta_map: Dict[int, float], plane: str) -> List[dict]:
    """ψ(s) ≈ Σ Δs / β(s) using β after each element (piecewise-constant)."""
    psi = 0.0
    out = []
    for el in elements:
        idx = el["ele_index"]
        beta = beta_map.get(idx)
        if beta is None or beta <= 0:
            # fall back to previous or skip advance
            beta_prev = beta_map.get(idx - 1, 1.0)
            beta = beta_prev if beta_prev > 0 else 1.0
        dpsi = el["L"] / beta
        psi_center = psi + 0.5 * dpsi
        e = dict(el)
        e[f"beta_{plane}"] = float(beta)
        e[f"psi_{plane}_rad"] = float(psi_center)
        e[f"psi_{plane}_turns"] = float(psi_center / (2 * np.pi))
        out.append(e)
        psi += dpsi
    return out


def attach_proxy(elements: List[dict], Qx: float, Qy: float) -> List[dict]:
    C = elements[-1]["s_end"] if elements else 1.0
    out = []
    for el in elements:
        frac = el["s_center"] / C if C else 0.0
        e = dict(el)
        e["mu_x_proxy_rad"] = 2 * np.pi * Qx * frac
        e["mu_y_proxy_rad"] = 2 * np.pi * Qy * frac
        out.append(e)
    return out


def wrap_pi(d: float) -> float:
    return float((d + np.pi) % (2 * np.pi) - np.pi)


def analyze_stem(stem: str) -> Optional[dict]:
    fox = STRUCTURES / stem / f"{stem}.fox"
    if not fox.is_file():
        return None
    text = fox.read_text(encoding="utf-8")
    lengths = parse_lengths(text)
    elements = element_index_map(walk_lattice(text, lengths))
    if not elements:
        return None

    bx = read_beta(DAT_SRC / stem / "BETAX")
    by = read_beta(DAT_SRC / stem / "BETAY")
    if not bx or not by:
        return {"stem": stem, "error": "missing BETAX/BETAY in COSY/src/dat"}

    df = get_data(stem, "FR0")
    if df is None:
        return {"stem": stem, "error": "missing FR0 mapping"}
    m = get_dual_model(df)
    Istar = solve_currents(m["R_int"], m["int_nat"], np.zeros(3))
    strengths = parse_sext_strengths_from_istar(Istar)
    Qx, Qy = default_tunes(stem)

    elements = integrate_phase(elements, bx, "x")
    # rebuild y on same list
    tmp = integrate_phase(elements, by, "y")
    for a, b in zip(elements, tmp):
        a["beta_y"] = b["beta_y"]
        a["psi_y_rad"] = b["psi_y_rad"]
        a["psi_y_turns"] = b["psi_y_turns"]
    elements = attach_proxy(elements, Qx, Qy)

    sexts = []
    diffs = []
    for el in elements:
        if el["kind"] != "MH" or el["family"] not in strengths:
            continue
        e = dict(el)
        e["K2"] = strengths[el["family"]]
        # fields expected by polar plot helpers
        e["mu_x_rad"] = e["psi_x_rad"]
        e["mu_y_rad"] = e["psi_y_rad"]
        dx = wrap_pi(e["psi_x_rad"] - e["mu_x_proxy_rad"])
        dy = wrap_pi(e["psi_y_rad"] - e["mu_y_proxy_rad"])
        e["dpsi_x_vs_proxy"] = dx
        e["dpsi_y_vs_proxy"] = dy
        sexts.append(e)
        diffs.append({"family": e["family"], "s": e["s_center"], "dx": dx, "dy": dy})

    res_x = resultant_vector(sexts, "x")
    res_y = resultant_vector(sexts, "y")
    res_x["note"] = "Σ |K2| e^{iψ}; ψ from Twiss ∫ds/β"
    res_y["note"] = "Σ |K2| e^{iψ}; ψ from Twiss ∫ds/β"

    # Temporary override of plot title via monkeypatch fields: plot_polar uses mu_* labels.
    # Write dedicated Twiss plots.
    px = _plot_twiss_polar(stem, sexts, "x", res_x)
    py = _plot_twiss_polar(stem, sexts, "y", res_y)

    rms_dx = float(np.sqrt(np.mean([d["dx"] ** 2 for d in diffs]))) if diffs else float("nan")
    rms_dy = float(np.sqrt(np.mean([d["dy"] ** 2 for d in diffs]))) if diffs else float("nan")

    out = {
        "stem": stem,
        "Qx_used": Qx,
        "Qy_used": Qy,
        "phase_source": "Twiss BETAX/BETAY from COSY/src/dat (linear optics; sextupoles do not change β on closed orbit)",
        "proxy_note": "proxy μ=2πQ·s/C for comparison only",
        "I_star": Istar.tolist(),
        "strengths": strengths,
        "n_sextupole_instances": len(sexts),
        "rms_dpsi_x_rad": rms_dx,
        "rms_dpsi_y_rad": rms_dy,
        "max_abs_dpsi_x_rad": float(max(abs(d["dx"]) for d in diffs)) if diffs else None,
        "max_abs_dpsi_y_rad": float(max(abs(d["dy"]) for d in diffs)) if diffs else None,
        "sextupoles": sexts,
        "resultant": {"x": res_x, "y": res_y},
        "plots": {"x": str(px), "y": str(py)},
    }
    write_json(DAT_OUT / stem / "sext_phase_twiss.json", out)
    return out


def _plot_twiss_polar(stem: str, sexts: List[dict], plane: str, resultant: dict) -> Path:
    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="polar")
    colors = {"SF1": "tab:red", "SF2": "tab:orange", "SD": "tab:blue"}
    radii = []
    for el in sexts:
        fam = el["family"]
        mu = el["mu_x_rad"] if plane == "x" else el["mu_y_rad"]
        ang = float(mu) + (np.pi if el["K2"] < 0 else 0.0)
        r = abs(float(el["K2"]))
        radii.append(r)
        ax.annotate(
            "",
            xy=(ang, r),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=colors.get(fam, "k"), lw=1.8),
        )
        ax.text(ang, r * 1.08, f"{fam}\ns={el['s_center']:.1f}", fontsize=7, ha="center")
    r_sum = float(resultant["magnitude"])
    a_sum = float(resultant["angle_rad"])
    radii.append(r_sum)
    ax.annotate(
        "",
        xy=(a_sum, r_sum),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color="black", lw=2.6, linestyle="--"),
    )
    rmax = max(radii) if radii else 1.0
    ax.set_rlim(0, rmax * 1.25 if rmax > 0 else 1.0)
    ax.set_title(
        f"{stem}: sextupole instances @ I*  (ψ_{plane}=∫ds/β)\n"
        f"|Σ|={resultant['magnitude']:.4g},  ∠Σ={resultant['angle_deg']:.1f}°",
        pad=18,
    )
    for fam, col in colors.items():
        ax.plot([], [], color=col, lw=2, label=fam)
    ax.plot([], [], color="black", lw=2.6, ls="--", label="resultant Σ")
    ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08), fontsize=8)
    out = PLOTS / f"{stem}_sext_phase_twiss_{plane}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    lines = []
    for stem in stems():
        out = analyze_stem(stem)
        if out is None or "error" in out:
            print(f"SKIP {stem}: {out}")
            continue
        lines.append(
            f"- {stem}: {out['n_sextupole_instances']} instances; "
            f"rms|Δψ|_x={out['rms_dpsi_x_rad']:.4g}, rms|Δψ|_y={out['rms_dpsi_y_rad']:.4g}; "
            f"|Σ|_x={out['resultant']['x']['magnitude']:.4g}, |Σ|_y={out['resultant']['y']['magnitude']:.4g}"
        )
        print(f"OK twiss-phase {stem}")

    # Keep geometric proxy script outputs, but document Twiss as primary.
    append_worklog(
        "## Twiss phase at sextupoles\n\n"
        "- **Статус:** verified (offline from existing BETAX/BETAY)\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_twiss_phase.py`\n"
        "- Linear β unchanged by sextupoles on closed orbit → existing Twiss usable at I*.\n"
        + "\n".join(lines)
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
