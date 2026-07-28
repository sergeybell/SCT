#!/usr/bin/env python3
"""Parse working-point COSY outputs and compute resonance distances."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DAT_OUT, PLOTS, append_worklog, load_config, stems, write_json  # noqa: E402

import matplotlib.pyplot as plt


def parse_wp_dat(path: Path) -> dict:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if k == "tag":
            out["tag"] = v
        else:
            out[k] = float(v.replace("D", "E"))
    return out


def resonance_distances(Qx: float, Qy: float, gammaG: float, P: int, ks: List[int]) -> List[dict]:
    rows = []
    for k in ks:
        for plane, Q in (("x", Qx), ("y", Qy)):
            for sign, s in (("plus", +1), ("minus", -1)):
                target = k * P + s * Q
                rows.append(
                    {
                        "plane": plane,
                        "k": k,
                        "form": f"{k}*P {'+' if s>0 else '-'} Q_{plane}",
                        "target": target,
                        "distance": abs(gammaG - target),
                    }
                )
    rows.sort(key=lambda r: r["distance"])
    return rows[:12]  # nearest


def offline_placeholder(stem: str, cfg: dict) -> dict:
    """If COSY not run yet, store expected structure from zero_point."""
    zp = json.loads((DAT_OUT / stem / "zero_point.json").read_text(encoding="utf-8"))
    P = int(cfg["superperiodicity"][stem])
    gG = float(cfg["gamma"]) * float(cfg["G"])
    points = {}
    for pt in zp["control_points"]["points"]:
        tag = pt["tag"]
        path = DAT_OUT / stem / f"wp_{tag}.dat"
        if path.is_file():
            points[tag] = parse_wp_dat(path)
            points[tag]["nearest_resonances"] = resonance_distances(
                points[tag]["Qx"], points[tag]["Qy"], points[tag]["gammaG"], P, cfg["resonance_orders_k"]
            )
        else:
            points[tag] = {
                "status": "pending_cosy",
                "I": pt["I"],
                "expected_gammaG": gG,
                "superperiodicity": P,
            }
    return {"stem": stem, "P": P, "gammaG_config": gG, "points": points}


def plot_wp_summary(all_data: Dict[str, dict]) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    for stem, data in all_data.items():
        for tag, pt in data["points"].items():
            if "Qx" not in pt:
                continue
            ax.scatter([pt["Qx"]], [pt["Qy"]], s=60, label=f"{stem}/{tag}")
            ax.annotate(f"{stem}:{tag}", (pt["Qx"], pt["Qy"]), fontsize=7)
    ax.set_xlabel("Qx")
    ax.set_ylabel("Qy")
    ax.set_title("Working points (from COSY WP runs)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    out = PLOTS / "working_points_summary.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    cfg = load_config()
    all_data = {}
    n_ready = 0
    for stem in stems():
        data = offline_placeholder(stem, cfg)
        write_json(DAT_OUT / stem / "working_point.json", data)
        all_data[stem] = data
        n_ready += sum(1 for p in data["points"].values() if "Qx" in p)
        print(f"OK working_point table {stem}: {sum(1 for p in data['points'].values() if 'Qx' in p)} COSY points")

    plot = None
    if n_ready:
        plot = plot_wp_summary(all_data)

    append_worklog(
        "## Working points\n\n"
        f"- **Статус:** {'partial' if n_ready else 'prepared'} ({n_ready} COSY wp dat parsed)\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_working_point.py`\n"
        "- После `run_cosy_jobs.py --job working_point` перезапустить этот скрипт.\n"
        f"- Plot: {plot.name if plot else 'n/a'}\n"
        f"- γG(config)={cfg['gamma']*cfg['G']:.6f}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
