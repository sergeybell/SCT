#!/usr/bin/env python3
"""Postprocess INJECT+TR outputs: per-particle Dbar, Delta nu_s, coherence C(n)."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DAT_OUT, PLOTS, append_worklog, load_config, write_json  # noqa: E402


def read_trpray(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                it = int(float(parts[0]))
                ray = int(float(parts[1]))
                vals = [float(x.replace("D", "E")) for x in parts[2:8]]
            except ValueError:
                continue
            rows.append([it, ray, *vals])
    if not rows:
        return pd.DataFrame(columns=["turn", "ray", "X", "A", "Y", "B", "T", "D"])
    return pd.DataFrame(rows, columns=["turn", "ray", "X", "A", "Y", "B", "T", "D"])


def read_trpspi(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                if len(parts) == 5:
                    turn, ray = int(float(parts[0])), int(float(parts[1]))
                    sx, sy, sz = [float(x.replace("D", "E")) for x in parts[2:5]]
                elif len(parts) >= 6:
                    turn = int(float(parts[0]))
                    ray = int(float(parts[2]))
                    sx, sy, sz = [float(x.replace("D", "E")) for x in parts[-3:]]
                else:
                    continue
            except ValueError:
                continue
            rows.append([turn, ray, sx, sy, sz])
    if not rows:
        return pd.DataFrame(columns=["turn", "ray", "Sx", "Sy", "Sz"])
    return pd.DataFrame(rows, columns=["turn", "ray", "Sx", "Sy", "Sz"])


def classify_rays(pray: Path) -> Dict[int, str]:
    """Map ray id -> group using initial PRAY if available."""
    groups = {}
    if not pray.is_file():
        return groups
    # PRAY dump format varies; fallback by injection order in analyze
    return groups


def assign_groups_from_initial(df: pd.DataFrame) -> Dict[int, str]:
    """At turn==min, classify by dominant coordinate."""
    t0 = df["turn"].min()
    init = df[df["turn"] == t0].sort_values("ray")
    groups = {}
    for _, row in init.iterrows():
        rid = int(row.ray)
        x, y, d = abs(row.X), abs(row.Y), abs(row.D)
        if x < 1e-12 and y < 1e-12 and d < 1e-12:
            groups[rid] = "ref"
        elif x >= y and x >= d:
            groups[rid] = "X"
        elif y >= x and y >= d:
            groups[rid] = "Y"
        else:
            groups[rid] = "D"
    return groups


def analyze_case(stem: str, tag: str, kind: str) -> Optional[dict]:
    base = DAT_OUT / stem
    trpray = base / f"track_{tag}_{kind}_TRPRAY.dat"
    trpspi = base / f"track_{tag}_{kind}_TRPSPI.dat"
    if not trpray.is_file() or not trpspi.is_file():
        return None

    orb = read_trpray(trpray)
    spi = read_trpspi(trpspi)
    if orb.empty or spi.empty:
        return {"error": "empty tracking files", "tag": tag, "kind": kind}

    groups = assign_groups_from_initial(orb)
    turns = np.sort(orb["turn"].unique())
    # drop first 20% as transient
    n1 = turns[int(0.2 * len(turns))] if len(turns) > 5 else turns[0]
    n2 = turns[-1]

    # Prefer reference with nonzero spin; skip COSY dummy ray 0 if |S|=0
    spin_norm = spi.groupby("ray")[["Sx", "Sy", "Sz"]].apply(
        lambda g: float(np.sqrt((g.iloc[0] ** 2).sum()))
    )
    ref_candidates = [r for r, g in groups.items() if g == "ref"]
    ref_id = None
    for r in ref_candidates:
        if spin_norm.get(r, 0.0) > 0.5:
            ref_id = r
            break
    if ref_id is None:
        # fallback: smallest ray with |S|~1 and near-zero orbit at t0
        t0 = orb["turn"].min()
        init = orb[orb["turn"] == t0]
        for _, row in init.sort_values("ray").iterrows():
            rid = int(row.ray)
            if spin_norm.get(rid, 0.0) > 0.5 and abs(row.X) + abs(row.Y) + abs(row.D) < 1e-12:
                ref_id = rid
                groups[rid] = "ref"
                break
    if ref_id is None:
        ref_id = int(orb["ray"].min())

    orb_ref = orb[orb["ray"] == ref_id]
    D_ref = float(orb_ref.loc[orb_ref["turn"] >= n1, "D"].mean())

    per_ray = []
    for rid, g in groups.items():
        sub = orb[orb["ray"] == rid]
        Dbar = float(sub.loc[sub["turn"] >= n1, "D"].mean())
        per_ray.append(
            {
                "ray": rid,
                "group": g,
                "Dbar": Dbar,
                "delta_eq": Dbar - D_ref,
                "X0": float(sub.loc[sub["turn"] == turns[0], "X"].iloc[0]),
                "Y0": float(sub.loc[sub["turn"] == turns[0], "Y"].iloc[0]),
                "D0": float(sub.loc[sub["turn"] == turns[0], "D"].iloc[0]),
            }
        )

    # spin phases in XZ plane (horizontal polarization diagnostic with Sy initial)
    # For PSI=90°, initial S=(0,1,0) — use phase in XZ relative to stable axis approx atan2(Sx,Sz)
    spi = spi.sort_values(["ray", "turn"])
    dnu = []
    phase_series = {}
    for rid in sorted(groups):
        s = spi[spi["ray"] == rid]
        if s.empty:
            continue
        phi = np.unwrap(np.arctan2(s["Sx"].values, s["Sz"].values))
        turns_s = s["turn"].values.astype(float)
        phase_series[rid] = (turns_s, phi)
        if len(turns_s) >= 2:
            # linear fit after transient
            mask = turns_s >= n1
            if mask.sum() >= 2:
                coeff = np.polyfit(turns_s[mask], phi[mask], 1)
                dnu.append({"ray": rid, "group": groups[rid], "dnu_s": coeff[0] / (2 * np.pi)})

    # coherence relative to ref
    if ref_id in phase_series:
        t_ref, phi_ref = phase_series[ref_id]
        # interpolate all to common turns
        common = turns
        dphi = []
        valid_rays = []
        for rid, (t, phi) in phase_series.items():
            if rid == ref_id:
                continue
            if groups.get(rid) == "ref":
                continue
            ip = np.interp(common, t, phi)
            ir = np.interp(common, t_ref, phi_ref)
            dphi.append(ip - ir)
            valid_rays.append(rid)
        dphi = np.array(dphi) if dphi else np.zeros((0, len(common)))
        if len(dphi):
            C = np.abs(np.mean(np.exp(1j * dphi), axis=0))
            rms_circ = np.sqrt(-2 * np.log(np.clip(C, 1e-15, 1.0)))
        else:
            C = np.ones(len(common))
            rms_circ = np.zeros(len(common))
    else:
        common, C, rms_circ = turns, np.ones(len(turns)), np.zeros(len(turns))

    summary = {
        "stem": stem,
        "tag": tag,
        "kind": kind,
        "n_rays": len(groups),
        "groups": groups,
        "D_ref": D_ref,
        "transient_cut_turn": float(n1),
        "per_ray": per_ray,
        "dnu_s": dnu,
        "delta_eq_rms": float(np.sqrt(np.mean([p["delta_eq"] ** 2 for p in per_ray if p["group"] != "ref"]))) if per_ray else None,
        "dnu_rms": float(np.sqrt(np.mean([d["dnu_s"] ** 2 for d in dnu]))) if dnu else None,
        "delta_eq_rms_by_group": {
            gname: float(np.sqrt(np.mean([p["delta_eq"] ** 2 for p in per_ray if p["group"] == gname] or [0.0])))
            for gname in ("X", "Y", "D")
        },
        "dnu_rms_by_group": {
            gname: float(np.sqrt(np.mean([d["dnu_s"] ** 2 for d in dnu if d["group"] == gname] or [0.0])))
            for gname in ("X", "Y", "D")
        },
        "C_final": float(C[-1]) if len(C) else None,
        "turns": common.tolist(),
        "C": C.tolist(),
        "rms_circular_phase": rms_circ.tolist(),
    }

    # plots
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle(f"{stem} / {tag} / {kind}")
    ax = axes[0, 0]
    for gname, color in [("X", "r"), ("Y", "b"), ("D", "g")]:
        xs = [p["X0"] ** 2 + p["Y0"] ** 2 + p["D0"] ** 2 for p in per_ray if p["group"] == gname]
        ys = [p["delta_eq"] for p in per_ray if p["group"] == gname]
        ax.scatter(xs, ys, c=color, label=gname)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"amp proxy $X_0^2+Y_0^2+D_0^2$")
    ax.set_ylabel(r"$\overline{D}_i-\overline{D}_{ref}$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    for d in dnu:
        ax.scatter([d["ray"]], [d["dnu_s"]], c={"X": "r", "Y": "b", "D": "g", "ref": "k"}.get(d["group"], "gray"))
    ax.set_xlabel("ray")
    ax.set_ylabel(r"$\Delta\nu_{s,i}$")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(common, C, lw=2)
    ax.set_xlabel("turn")
    ax.set_ylabel("C(n)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(common, rms_circ, lw=2, color="purple")
    ax.set_xlabel("turn")
    ax.set_ylabel("RMS circular phase")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pout = PLOTS / f"{stem}_track_{tag}_{kind}.png"
    fig.savefig(pout, dpi=150, bbox_inches="tight")
    plt.close(fig)
    summary["plot"] = str(pout)
    write_json(base / f"track_{tag}_{kind}_analysis.json", summary)
    return summary


def compare_plot(stem: str, kind: str, tags: list) -> Optional[Path]:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(f"{stem}: tracking compare ({kind})")
    any_data = False
    xg, yg = [], []
    for tag in tags:
        path = DAT_OUT / stem / f"track_{tag}_{kind}_analysis.json"
        if not path.is_file():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        any_data = True
        axes[0].plot(data["turns"], data["C"], lw=2, label=tag)
        by = data.get("delta_eq_rms_by_group", {})
        xg.append(by.get("X", 0.0))
        yg.append(by.get("Y", 0.0))
        axes[1].bar(tag, by.get("X", 0.0))
        axes[2].bar(tag, by.get("Y", 0.0))
    if not any_data:
        plt.close(fig)
        return None
    axes[0].set_xlabel("turn")
    axes[0].set_ylabel("C(n)")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    axes[1].set_ylabel(r"RMS $\Delta\delta_{eq}$ (X group)")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].grid(True, alpha=0.3)
    axes[2].set_ylabel(r"RMS $\Delta\delta_{eq}$ (Y group)")
    axes[2].tick_params(axis="x", rotation=30)
    axes[2].grid(True, alpha=0.3)
    out = PLOTS / f"{stem}_track_compare_{kind}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    cfg = load_config()
    stem = cfg["tracking"]["pilot_stem"]
    tags = ["natural", "Istar", "ctrl_xi_x", "ctrl_xi_y", "ctrl_eta1"]
    done = []
    for kind in ("smoke", "full"):
        for tag in tags:
            s = analyze_case(stem, tag, kind)
            if s and "error" not in s:
                done.append(f"{tag}/{kind}: C_final={s.get('C_final')}, dde_rms={s.get('delta_eq_rms')}")
                print(f"OK analyze {tag} {kind}")
            elif s is None:
                print(f"SKIP missing {tag} {kind}")
        compare_plot(stem, kind, tags)

    append_worklog(
        "## Tracking analysis\n\n"
        f"- **Статус:** {'verified' if done else 'prepared (waiting COSY track outputs)'}\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_tracking.py`\n"
        + ("\n".join(f"- {x}" for x in done) if done else "- No TRPRAY/TRPSPI yet.\n")
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
