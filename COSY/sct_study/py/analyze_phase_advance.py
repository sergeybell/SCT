#!/usr/bin/env python3
"""Parse lattice FOX → physical sextupole s-positions; plot 2D polar phase diagrams."""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DAT_OUT,
    PLOTS,
    STRUCTURES,
    append_worklog,
    get_data,
    get_dual_model,
    load_config,
    solve_currents,
    stems,
    write_json,
)


def parse_lengths(text: str) -> Dict[str, float]:
    vals = {}
    for m in re.finditer(r"(L_\w+)\s*:=\s*([0-9.eE+-]+)", text):
        vals[m.group(1)] = float(m.group(2))
    return vals


def parse_sext_strengths_from_istar(Istar: np.ndarray) -> Dict[str, float]:
    # SF1<-SGx1, SF2<-SGx2, SD<-SGy1
    return {"SF1": float(Istar[0]), "SF2": float(Istar[1]), "SD": float(Istar[2])}


def walk_lattice(text: str, lengths: Dict[str, float]) -> List[dict]:
    body = text.split("{BEGIN LATTICE}")[-1].split("ENDLOOP")[0]
    s = 0.0
    mux = 0.0  # placeholder: geometric phase proxy via s/L * 2πQx later
    elements = []
    # element patterns
    patterns = [
        (r"DL\s+(L_\w+)", "DL"),
        (r"QUAD\s+(L_\w+)", "QUAD"),
        (r"SBEND\s+(L_\w+)", "SBEND"),
        (r"MH\s+(L_\w+)\s+(\w+)", "MH"),
    ]
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith("{"):
            continue
        for pat, kind in patterns:
            m = re.search(pat, line)
            if not m:
                continue
            Lname = m.group(1)
            L = lengths.get(Lname, 0.0)
            fam = None
            if kind == "MH":
                fam = m.group(2)  # SF1/SF2/SD
                comment = re.search(r"\{(\w+)\}", line)
                if comment:
                    fam = comment.group(1)
            elements.append(
                {
                    "kind": kind,
                    "L": L,
                    "s_center": s + 0.5 * L,
                    "s_end": s + L,
                    "family": fam,
                    "Lname": Lname,
                }
            )
            s += L
            break
    return elements


def estimate_phase_advances(elements: List[dict], Qx: float, Qy: float) -> List[dict]:
    """Uniform phase advance proxy μ = 2π Q * (s/C). Replace later with Twiss BETS if available."""
    C = elements[-1]["s_end"] if elements else 1.0
    out = []
    for el in elements:
        frac = el["s_center"] / C if C else 0.0
        e = dict(el)
        e["mu_x_rad"] = 2 * np.pi * Qx * frac
        e["mu_y_rad"] = 2 * np.pi * Qy * frac
        e["mu_x_turns"] = Qx * frac
        e["mu_y_turns"] = Qy * frac
        out.append(e)
    return out


def default_tunes(stem: str) -> Tuple[float, float]:
    # Placeholder fractional tunes until COSY WP run fills real values.
    # Use stem-dependent mild defaults; overwritten from working_point.json if present.
    wp = DAT_OUT / stem / "working_point.json"
    if wp.is_file():
        import json

        data = json.loads(wp.read_text(encoding="utf-8"))
        for tag in ("Istar", "natural"):
            if tag in data.get("points", {}):
                return float(data["points"][tag]["Qx"]), float(data["points"][tag]["Qy"])
    # OptiM-like rough guesses from period count (not physics-accurate)
    return 3.6 + 0.1 * int(stem.split("_")[1]), 3.2 + 0.05 * int(stem.split("_")[1])


def signed_polar(K: float, mu_rad: float) -> Tuple[float, float]:
    """Return (angle, radius) with negative K₂ rotated by π."""
    ang = float(mu_rad) + (np.pi if K < 0 else 0.0)
    return ang, abs(float(K))


def resultant_vector(sexts: List[dict], plane: str) -> dict:
    """Complex sum of signed sextupole arrows in the chosen plane."""
    z = 0.0 + 0.0j
    for el in sexts:
        mu = el["mu_x_rad"] if plane == "x" else el["mu_y_rad"]
        ang, r = signed_polar(el["K2"], mu)
        z += r * np.exp(1j * ang)
    return {
        "plane": plane,
        "Re": float(np.real(z)),
        "Im": float(np.imag(z)),
        "magnitude": float(np.abs(z)),
        "angle_rad": float(np.angle(z)),
        "angle_deg": float(np.degrees(np.angle(z))),
        "note": "Diagnostic sum Σ |K2| e^{iμ_signed}; μ is s/C·Q proxy, not Twiss phase",
    }


def plot_polar(stem: str, sexts: List[dict], plane: str, resultant: dict) -> Path:
    """One figure per stem/plane with its own radial scale and resultant arrow."""
    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="polar")
    colors = {"SF1": "tab:red", "SF2": "tab:orange", "SD": "tab:blue"}

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
            arrowprops=dict(arrowstyle="->", color=colors.get(fam, "k"), lw=1.8),
        )
        ax.text(ang, r * 1.08, f"{fam}\ns={el['s_center']:.1f} m", fontsize=7, ha="center", va="bottom")

    # Resultant vector (contrasting style)
    r_sum = float(resultant["magnitude"])
    a_sum = float(resultant["angle_rad"])
    radii.append(r_sum)
    ax.annotate(
        "",
        xy=(a_sum, r_sum),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color="black", lw=2.6, linestyle="--"),
    )
    ax.text(
        a_sum,
        max(r_sum, 1e-12) * 1.12,
        "Σ",
        fontsize=11,
        fontweight="bold",
        ha="center",
        color="black",
    )

    # Own radial scale for this structure (gradients differ between stems)
    rmax = max(radii) if radii else 1.0
    ax.set_rlim(0, rmax * 1.25 if rmax > 0 else 1.0)
    ax.set_title(
        f"{stem}: sextupoles @ I*  (μ_{plane}=2πQ·s/C proxy)\n"
        f"|Σ|={resultant['magnitude']:.4g},  ∠Σ={resultant['angle_deg']:.1f}°",
        pad=18,
    )
    # Legend proxies
    for fam, col in colors.items():
        ax.plot([], [], color=col, lw=2, label=fam)
    ax.plot([], [], color="black", lw=2.6, ls="--", label="resultant Σ")
    ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08), fontsize=8)

    out = PLOTS / f"{stem}_sext_phase_{plane}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    lines = []
    for stem in stems():
        text = (STRUCTURES / stem / f"{stem}.fox").read_text(encoding="utf-8")
        lengths = parse_lengths(text)
        elements = walk_lattice(text, lengths)
        df = get_data(stem, "FR0")
        if df is None:
            continue
        m = get_dual_model(df)
        Istar = solve_currents(m["R_int"], m["int_nat"], np.zeros(3))
        strengths = parse_sext_strengths_from_istar(Istar)
        Qx, Qy = default_tunes(stem)
        elements = estimate_phase_advances(elements, Qx, Qy)
        sexts = []
        for el in elements:
            if el["kind"] != "MH" or el["family"] not in strengths:
                continue
            e = dict(el)
            e["K2"] = strengths[el["family"]]
            sexts.append(e)
        res_x = resultant_vector(sexts, "x")
        res_y = resultant_vector(sexts, "y")
        write_json(
            DAT_OUT / stem / "sext_phase.json",
            {
                "stem": stem,
                "Qx_used": Qx,
                "Qy_used": Qy,
                "phase_note": "Uniform μ=2πQ·s/C proxy until Twiss-based phases available",
                "I_star": Istar.tolist(),
                "strengths": strengths,
                "sextupoles": sexts,
                "resultant": {"x": res_x, "y": res_y},
            },
        )
        px = plot_polar(stem, sexts, "x", res_x)
        py = plot_polar(stem, sexts, "y", res_y)
        lines.append(
            f"- {stem}: {len(sexts)} magnets; |Σ|_x={res_x['magnitude']:.4g}, "
            f"|Σ|_y={res_y['magnitude']:.4g}; {px.name}, {py.name}"
        )
        print(f"OK phase {stem}: {len(sexts)} sextupoles")

    append_worklog(
        "## Phase diagrams\n\n"
        "- **Статус:** prepared (s/C·Q proxy for μ; update after working_point COSY run)\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_phase_advance.py`\n"
        "- Диаграмма: длина = |K₂|, угол = μ; знак K₂ → +π; чёрный пунктир = результирующий Σ.\n"
        "- Масштаб радиуса — свой для каждой структуры.\n"
        + "\n".join(lines) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
