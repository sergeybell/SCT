#!/usr/bin/env python3
"""Postprocess INJECT+TR: mean_D_offset, spin coherence, relative spin tunes."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


def assign_groups_from_inject(orb: pd.DataFrame, num_per_group: int) -> Tuple[Dict[int, str], int]:
    """Classify by INJECT layout: ray0 dummy, ray1 ref, then X/Y/D blocks.

    Cosy always inserts ray 0. Our template then injects one reference plus
    three groups of size ``num_per_group``. Classification by |X|,|Y|,|D| is
    wrong when linspace includes a zero amplitude ray.
    """
    rays = sorted(int(r) for r in orb["ray"].unique())
    groups: Dict[int, str] = {}
    if not rays:
        return groups, 0

    # Prefer explicit layout starting at ray 0.
    ref_id = 1 if 1 in rays else rays[0]
    groups[0] = "dummy" if 0 in rays else groups.get(0, "dummy")
    if 0 in rays:
        groups[0] = "dummy"
    groups[ref_id] = "ref"

    x0 = ref_id + 1
    for i in range(num_per_group):
        groups[x0 + i] = "X"
    y0 = x0 + num_per_group
    for i in range(num_per_group):
        groups[y0 + i] = "Y"
    d0 = y0 + num_per_group
    for i in range(num_per_group):
        groups[d0 + i] = "D"

    # Any unexpected ray ids: mark unknown (still tracked).
    for r in rays:
        groups.setdefault(r, "other")
    return groups, ref_id


def nyquist_turns(nu_s: float) -> float:
    """Max save interval (turns) for unambiguous spin-tune sampling."""
    return 0.5 / max(abs(nu_s), 1e-12)


def median_save_interval(turns: np.ndarray) -> float:
    if len(turns) < 2:
        return float("nan")
    return float(np.median(np.diff(np.sort(turns.astype(float)))))


def fit_phase_slope(
    turns: np.ndarray, phi: np.ndarray, n1: float
) -> Tuple[Optional[float], Optional[float], int]:
    mask = turns >= n1
    if mask.sum() < 2:
        return None, None, 0
    t = turns[mask]
    p = phi[mask]
    coeff = np.polyfit(t, p, 1)
    resid = p - np.polyval(coeff, t)
    rms = float(np.sqrt(np.mean(resid**2)))
    return float(coeff[0] / (2 * np.pi)), rms, int(mask.sum())


def analyze_case(stem: str, tag: str, kind: str, cfg: Optional[dict] = None) -> Optional[dict]:
    cfg = cfg or load_config()
    tr_cfg = cfg["tracking"]
    num = int(tr_cfg["num_per_group"])
    nu_s_ref = abs(float(cfg["gamma"]) * float(cfg["G"]))

    base = DAT_OUT / stem
    trpray = base / f"track_{tag}_{kind}_TRPRAY.dat"
    trpspi = base / f"track_{tag}_{kind}_TRPSPI.dat"
    if not trpray.is_file() or not trpspi.is_file():
        return None

    orb = read_trpray(trpray)
    spi = read_trpspi(trpspi)
    if orb.empty or spi.empty:
        return {"error": "empty tracking files", "tag": tag, "kind": kind}

    groups, ref_id = assign_groups_from_inject(orb, num)
    turns = np.sort(orb["turn"].unique())
    n1 = turns[int(0.2 * len(turns))] if len(turns) > 5 else turns[0]
    save_dt = median_save_interval(turns)
    nyq = nyquist_turns(nu_s_ref)
    spin_aliased = bool(save_dt > nyq)
    spin_samples = int(len(turns))
    spin_turns_physical = float(nu_s_ref * (turns[-1] - turns[0])) if len(turns) else 0.0

    orb_ref = orb[orb["ray"] == ref_id]
    D_ref = float(orb_ref.loc[orb_ref["turn"] >= n1, "D"].mean()) if not orb_ref.empty else float("nan")

    per_ray = []
    for rid, g in groups.items():
        if g in ("dummy", "other"):
            continue
        sub = orb[orb["ray"] == rid]
        if sub.empty:
            continue
        Dbar = float(sub.loc[sub["turn"] >= n1, "D"].mean())
        per_ray.append(
            {
                "ray": rid,
                "group": g,
                "Dbar": Dbar,
                "mean_D_offset": Dbar - D_ref,
                "X0": float(sub.loc[sub["turn"] == turns[0], "X"].iloc[0]),
                "Y0": float(sub.loc[sub["turn"] == turns[0], "Y"].iloc[0]),
                "D0": float(sub.loc[sub["turn"] == turns[0], "D"].iloc[0]),
            }
        )

    spi = spi.sort_values(["ray", "turn"])
    phase_series: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    spin_rows = []
    for rid, g in groups.items():
        if g in ("dummy", "other"):
            continue
        s = spi[spi["ray"] == rid]
        if s.empty:
            continue
        sx, sy, sz = s["Sx"].values, s["Sy"].values, s["Sz"].values
        phi = np.unwrap(np.arctan2(sx, sz))
        turns_s = s["turn"].values.astype(float)
        phase_series[rid] = (turns_s, phi)
        slope, rms, nfit = fit_phase_slope(turns_s, phi, float(n1))
        pol = float(np.mean(np.sqrt(sx**2 + sy**2 + sz**2)))
        row = {
            "ray": rid,
            "group": g,
            "mean_polarization": pol,
            "phase_fit_rms_rad": rms,
            "n_fit": nfit,
        }
        if spin_aliased:
            row["aliased_phase_slope"] = slope
            row["dnu_s_rel"] = None
        else:
            row["aliased_phase_slope"] = None
            row["dnu_s_rel"] = slope  # absolute slope; relative filled below
        spin_rows.append(row)

    # Relative spin tune / phase vs reference
    if ref_id in phase_series:
        t_ref, phi_ref = phase_series[ref_id]
        common = turns
        dphi = []
        valid_rays = []
        for rid, (t, phi) in phase_series.items():
            if rid == ref_id or groups.get(rid) == "ref":
                continue
            ip = np.interp(common, t, phi)
            ir = np.interp(common, t_ref, phi_ref)
            dphi.append(ip - ir)
            valid_rays.append(rid)
            # relative slope after transient
            rel_slope, rel_rms, nfit = fit_phase_slope(
                common.astype(float), ip - ir, float(n1)
            )
            for row in spin_rows:
                if row["ray"] == rid:
                    if spin_aliased:
                        row["aliased_phase_slope_rel"] = rel_slope
                        row["dnu_s_rel"] = None
                    else:
                        row["dnu_s_rel"] = rel_slope
                        row["phase_fit_rms_rad"] = rel_rms
                    break
        dphi_arr = np.array(dphi) if dphi else np.zeros((0, len(common)))
        if len(dphi_arr):
            C = np.abs(np.mean(np.exp(1j * dphi_arr), axis=0))
            rms_circ = np.sqrt(-2 * np.log(np.clip(C, 1e-15, 1.0)))
        else:
            C = np.ones(len(common))
            rms_circ = np.zeros(len(common))
        # mean polarization of ensemble (excluding dummy)
        ens = spi[spi["ray"].isin([r for r, g in groups.items() if g not in ("dummy", "other")])]
        if not ens.empty:
            mag = np.sqrt(ens["Sx"] ** 2 + ens["Sy"] ** 2 + ens["Sz"] ** 2)
            # average vector polarization over particles at each saved turn
            pol_t = []
            for tt in common:
                chunk = ens[ens["turn"] == tt]
                if chunk.empty:
                    pol_t.append(float("nan"))
                    continue
                mx = chunk["Sx"].mean()
                my = chunk["Sy"].mean()
                mz = chunk["Sz"].mean()
                pol_t.append(float(np.sqrt(mx**2 + my**2 + mz**2)))
        else:
            pol_t = [float("nan")] * len(common)
    else:
        common, C, rms_circ = turns, np.ones(len(turns)), np.zeros(len(turns))
        pol_t = [float("nan")] * len(turns)

    def rms_of(key: str, group: Optional[str] = None) -> Optional[float]:
        vals = [
            p[key]
            for p in per_ray
            if p["group"] != "ref" and (group is None or p["group"] == group)
        ]
        return float(np.sqrt(np.mean(np.square(vals)))) if vals else None

    def rms_spin(key: str, group: Optional[str] = None) -> Optional[float]:
        vals = []
        for d in spin_rows:
            if d["group"] == "ref":
                continue
            if group is not None and d["group"] != group:
                continue
            v = d.get(key)
            if v is not None:
                vals.append(v)
        return float(np.sqrt(np.mean(np.square(vals)))) if vals else None

    summary = {
        "stem": stem,
        "tag": tag,
        "kind": kind,
        "n_rays": len([g for g in groups.values() if g not in ("dummy", "other")]),
        "groups": {str(k): v for k, v in groups.items()},
        "ref_ray": ref_id,
        "D_ref": D_ref,
        "transient_cut_turn": float(n1),
        "sampling": {
            "save_interval_turns": save_dt,
            "nyquist_max_interval": nyq,
            "nu_s_ref": nu_s_ref,
            "spin_aliased": spin_aliased,
            "n_saved_turns": spin_samples,
            "physical_spin_turns": spin_turns_physical,
            "observable_spin_turns_note": (
                "Physical spin revolutions ν_s·ΔN exist, but phase samples only "
                "resolve them if save_interval ≤ 1/(2|ν_s|)."
            ),
        },
        "nomenclature": {
            "D": "COSY TRPRAY 6th coordinate (relative momentum/energy deviation)",
            "mean_D_offset": "⟨D⟩_i − ⟨D⟩_ref after transient cut; tracking proxy, not theoretical Δδ_eq",
            "delta_eq_theory": "Senichev equilibrium momentum shift from ξ, η1, emittances",
            "C_n": "Phase coherence |⟨exp(iΔφ)⟩| at saved turns only",
            "dnu_s_rel": "Relative spin-tune estimate from unwrapped phase; None if aliased",
            "aliased_phase_slope": "Raw phase slope / 2π when sampling is aliased — not a physical tune",
            "mean_polarization": "|⟨S⟩| ensemble polarization (distinct from phase coherence)",
            "spin_plane": (
                "Phase = unwrap(atan2(Sx,Sz)); requires horizontal/longitudinal initial spin "
                "(psi_deg≈0). Vertical spin (psi=90°) makes this phase undefined."
            ),
        },
        "per_ray": per_ray,
        "spin_per_ray": spin_rows,
        "mean_D_offset_rms": rms_of("mean_D_offset"),
        "mean_D_offset_rms_by_group": {
            gname: rms_of("mean_D_offset", gname) for gname in ("X", "Y", "D")
        },
        "dnu_rms": rms_spin("dnu_s_rel") if not spin_aliased else None,
        "aliased_slope_rms": rms_spin("aliased_phase_slope") if spin_aliased else None,
        "dnu_rms_by_group": {
            gname: (rms_spin("dnu_s_rel") if not spin_aliased else rms_spin("aliased_phase_slope_rel"))
            for gname in ("X", "Y", "D")
        },
        "C_final": float(C[-1]) if len(C) else None,
        "polarization_final": float(pol_t[-1]) if pol_t else None,
        "turns": common.tolist(),
        "C": C.tolist(),
        "rms_circular_phase": rms_circ.tolist(),
        "mean_polarization_vs_turn": pol_t,
    }

    # plots
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    title = f"{stem} / {tag} / {kind}"
    if spin_aliased:
        title += "  [spin sampling ALIASED]"
    fig.suptitle(title)

    ax = axes[0, 0]
    for gname, color in [("X", "r"), ("Y", "b"), ("D", "g")]:
        xs = [p["X0"] ** 2 + p["Y0"] ** 2 + p["D0"] ** 2 for p in per_ray if p["group"] == gname]
        ys = [p["mean_D_offset"] for p in per_ray if p["group"] == gname]
        ax.scatter(xs, ys, c=color, label=gname)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"amp proxy $X_0^2+Y_0^2+D_0^2$")
    ax.set_ylabel(r"$\overline{D}_i-\overline{D}_{ref}$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    key = "aliased_phase_slope" if spin_aliased else "dnu_s_rel"
    for d in spin_rows:
        v = d.get(key)
        if v is None and spin_aliased:
            v = d.get("aliased_phase_slope_rel")
        if v is None:
            continue
        ax.scatter(
            [d["ray"]],
            [v],
            c={"X": "r", "Y": "b", "D": "g", "ref": "k"}.get(d["group"], "gray"),
        )
    ax.set_xlabel("ray")
    ax.set_ylabel("aliased slope" if spin_aliased else r"$\Delta\nu_{s,i}$ (rel)")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(common, C, lw=2, label="C(n)")
    ax.plot(common, pol_t, lw=1.5, ls="--", label=r"$|\langle S\rangle|$")
    ax.set_xlabel("turn")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8)
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
    for tag in tags:
        path = DAT_OUT / stem / f"track_{tag}_{kind}_analysis.json"
        if not path.is_file():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        any_data = True
        axes[0].plot(data["turns"], data["C"], lw=2, label=tag)
        by = data.get("mean_D_offset_rms_by_group") or data.get("delta_eq_rms_by_group", {})
        axes[1].bar(tag, by.get("X", 0.0) or 0.0)
        axes[2].bar(tag, by.get("Y", 0.0) or 0.0)
    if not any_data:
        plt.close(fig)
        return None
    axes[0].set_xlabel("turn")
    axes[0].set_ylabel("C(n)")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    axes[1].set_ylabel(r"RMS mean_D_offset (X)")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].grid(True, alpha=0.3)
    axes[2].set_ylabel(r"RMS mean_D_offset (Y)")
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
    kinds = ["smoke", "full", "dense"]
    done = []
    for kind in kinds:
        for tag in tags:
            s = analyze_case(stem, tag, kind, cfg)
            if s and "error" not in s:
                samp = s["sampling"]
                done.append(
                    f"{tag}/{kind}: C_final={s.get('C_final')}, "
                    f"mean_D_rms={s.get('mean_D_offset_rms')}, "
                    f"aliased={samp['spin_aliased']}, Δn={samp['save_interval_turns']}"
                )
                print(f"OK analyze {tag} {kind} aliased={samp['spin_aliased']}")
            elif s is None:
                print(f"SKIP missing {tag} {kind}")
        compare_plot(stem, kind, tags)

    append_worklog(
        "## Tracking analysis\n\n"
        f"- **Статус:** {'verified' if done else 'prepared (waiting COSY track outputs)'}\n"
        "- **Команда:** `python COSY/sct_study/py/analyze_tracking.py`\n"
        "- Метрики: `mean_D_offset` (не Δδ_eq), `C(n)`, `|⟨S⟩|`, "
        "`dnu_s_rel` только при save_interval ≤ 1/(2|ν_s|).\n"
        + ("\n".join(f"- {x}" for x in done) if done else "- No TRPRAY/TRPSPI yet.\n")
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
