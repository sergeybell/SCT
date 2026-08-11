#!/usr/bin/env python3
"""Spin resonance distance scan: imperfection, intrinsic, synchrotron sidebands."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DAT_OUT, PLOTS, append_worklog, load_config, stems, write_json  # noqa: E402


def frac(x: float) -> float:
    """Distance to nearest integer."""
    return abs(x - round(x))


def resonance_catalog(
    Qx: float,
    Qy: float,
    nu_s: float,
    Qs: float,
    *,
    ks: List[int],
    ms: List[int],
) -> List[dict]:
    """Build geometric resonance list for deuteron (and comparable proton) scan.

    Lines:
      imperfection:  ν_s = k
      intrinsic:     ν_s = k ± Q_x, k ± Q_y
      synchrotron:   ν_s = k ± m Q_s
      combined:      ν_s = k ± Q_{x,y} ± m Q_s
    Distance = min_n |ν_s − target − n| (fractional).
    """
    rows: List[dict] = []
    nu = float(nu_s)

    def add(kind: str, form: str, target: float, **extra):
        d = frac(nu - target)
        rows.append(
            {
                "kind": kind,
                "form": form,
                "target": float(target),
                "distance": float(d),
                **extra,
            }
        )

    for k in ks:
        add("imperfection", f"ν_s={k}", float(k), k=k, m=0, plane=None)
        for m in ms:
            if m == 0:
                continue
            for s, sm in ((+1, "+"), (-1, "-")):
                add(
                    "synchrotron",
                    f"ν_s={k}{sm}{m}Q_s",
                    k + s * m * Qs,
                    k=k,
                    m=m,
                    plane="s",
                )
        for plane, Q in (("x", Qx), ("y", Qy)):
            for s, sm in ((+1, "+"), (-1, "-")):
                add(
                    "intrinsic",
                    f"ν_s={k}{sm}Q_{plane}",
                    k + s * Q,
                    k=k,
                    m=0,
                    plane=plane,
                )
                for m in ms:
                    if m == 0:
                        continue
                    for s2, sm2 in ((+1, "+"), (-1, "-")):
                        add(
                            "combined",
                            f"ν_s={k}{sm}Q_{plane}{sm2}{m}Q_s",
                            k + s * Q + s2 * m * Qs,
                            k=k,
                            m=m,
                            plane=plane,
                        )
    rows.sort(key=lambda r: r["distance"])
    return rows


def proton_comparison(gamma: float) -> dict:
    """Reference numbers: same lattice energy would give larger |Gγ| for protons."""
    G_p = 1.7928474
    G_d = -0.142987
    return {
        "deuteron_G": G_d,
        "proton_G": G_p,
        "gamma": gamma,
        "deuteron_gammaG": gamma * G_d,
        "proton_gammaG_same_gamma": gamma * G_p,
        "note": (
            "Deuteron |Gγ|≪1 keeps low-order imperfection lines farther in absolute "
            "spin tune than for protons at the same γ; this does not prove resonances "
            "are irrelevant for coherence — check distances and observed excitation."
        ),
    }


def load_wp_point(stem: str, tag: str) -> Optional[dict]:
    path = DAT_OUT / stem / "working_point.json"
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    pt = data.get("points", {}).get(tag)
    if not pt or "Qx" not in pt:
        return None
    return pt


def analyze_stem(stem: str, cfg: dict) -> dict:
    ks = list(cfg["resonance_orders_k"])
    ms = list(cfg.get("synchrotron_sideband_m", [0, 1, 2]))
    Qs = float(cfg.get("Qs_estimate", {}).get("default", 0.01))
    tags = ["natural", "Istar", "ctrl_xi_x", "ctrl_xi_y", "ctrl_eta1"]
    points = {}
    for tag in tags:
        pt = load_wp_point(stem, tag)
        if pt is None:
            points[tag] = {"status": "missing_wp"}
            continue
        nu = float(pt.get("NU_s", abs(pt.get("gammaG", cfg["gamma"] * cfg["G"]))))
        # Use signed gammaG magnitude consistently with working-point NU_s (>0 in dat)
        catalog = resonance_catalog(
            float(pt["Qx"]),
            float(pt["Qy"]),
            nu,
            Qs,
            ks=ks,
            ms=ms,
        )
        nearest = {}
        for kind in ("imperfection", "intrinsic", "synchrotron", "combined"):
            sub = [r for r in catalog if r["kind"] == kind]
            nearest[kind] = sub[0] if sub else None
        points[tag] = {
            "Qx": pt["Qx"],
            "Qy": pt["Qy"],
            "NU_s": nu,
            "gammaG": pt.get("gammaG"),
            "Qs_used": Qs,
            "nearest_by_kind": nearest,
            "nearest_all": catalog[:15],
            "n_lines_scanned": len(catalog),
        }
    return {
        "stem": stem,
        "Qs_note": cfg.get("Qs_estimate", {}).get("note"),
        "proton_comparison": proton_comparison(float(cfg["gamma"])),
        "points": points,
    }


def plot_summary(all_data: Dict[str, dict]) -> Path:
    rows = []
    for stem, data in all_data.items():
        for tag, pt in data["points"].items():
            if "nearest_by_kind" not in pt:
                continue
            for kind, r in pt["nearest_by_kind"].items():
                if not r:
                    continue
                rows.append(
                    {
                        "stem": stem,
                        "tag": tag,
                        "kind": kind,
                        "distance": r["distance"],
                        "form": r["form"],
                    }
                )
    fig, ax = plt.subplots(figsize=(10, 5))
    if rows:
        df = pd.DataFrame(rows)
        # focus on Istar
        sub = df[df["tag"] == "Istar"]
        for kind, marker in [
            ("imperfection", "o"),
            ("intrinsic", "s"),
            ("synchrotron", "^"),
            ("combined", "D"),
        ]:
            ss = sub[sub["kind"] == kind]
            if ss.empty:
                continue
            ax.scatter(ss["stem"], ss["distance"], marker=marker, s=70, label=kind)
        ax.set_ylabel("fractional distance to nearest line")
        ax.set_title("Nearest spin resonances @ I* (geometric)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    out = PLOTS / "resonance_scan_summary.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def strength_proxy_from_tracking(stem: str, tag: str, kind: str = "dense") -> Optional[dict]:
    """If dense tracking analysis exists, report observed spin metrics as excitation proxy."""
    path = DAT_OUT / stem / f"track_{tag}_{kind}_analysis.json"
    if not path.is_file():
        # fall back to full/smoke
        for k in ("full", "smoke"):
            path = DAT_OUT / stem / f"track_{tag}_{k}_analysis.json"
            if path.is_file():
                kind = k
                break
        else:
            return None
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        "track_kind": kind,
        "C_final": data.get("C_final"),
        "polarization_final": data.get("polarization_final"),
        "mean_D_offset_rms_by_group": data.get("mean_D_offset_rms_by_group"),
        "sampling_aliased": data.get("sampling", {}).get("spin_aliased"),
        "dnu_rms": data.get("dnu_rms"),
        "note": (
            "Tracking metrics are observable proxies for decoherence/excitation; "
            "not Froissart–Stora resonance strength."
        ),
    }


def main() -> int:
    cfg = load_config()
    all_data = {}
    lines = []
    for stem in stems():
        data = analyze_stem(stem, cfg)
        # attach tracking proxies for pilot tags when available
        for tag, pt in data["points"].items():
            if isinstance(pt, dict) and "Qx" in pt:
                pt["tracking_proxy"] = strength_proxy_from_tracking(stem, tag)
        write_json(DAT_OUT / stem / "resonance_scan.json", data)
        all_data[stem] = data
        ist = data["points"].get("Istar", {})
        if "nearest_by_kind" in ist:
            n = ist["nearest_by_kind"]
            lines.append(
                f"- {stem}@Istar: imp={n['imperfection']['distance']:.4g} ({n['imperfection']['form']}), "
                f"intr={n['intrinsic']['distance']:.4g} ({n['intrinsic']['form']}), "
                f"sync={n['synchrotron']['distance']:.4g}, comb={n['combined']['distance']:.4g}"
            )
        print(f"OK resonance scan {stem}")

    write_json(DAT_OUT / "resonance_scan_summary.json", {"stems": all_data})
    plot = plot_summary(all_data)
    append_worklog(
        "## Resonance scan\n\n"
        "- **Статус:** verified (geometric distances + tracking proxies when present)\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_resonances.py`\n"
        f"- Plot: `{plot.name}`\n"
        "- Qs from config estimate until RF one-turn extraction is added.\n"
        + "\n".join(lines)
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
