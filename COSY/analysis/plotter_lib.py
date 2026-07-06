"""Twiss plotting utilities (extracted from Plotter.ipynb / twiss_plot.ipynb)."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd

FILES_TO_PLOT: Dict[str, Dict[str, str]] = {
    "BETAX": {"ax": "left", "color": "red", "label": r"$\beta_x$"},
    "BETAY": {"ax": "left", "color": "green", "label": r"$\beta_y$"},
    "DISPX": {"ax": "right", "color": "blue", "label": r"$D_x$"},
    "DISPY": {"ax": "right", "color": "black", "label": r"$D_y$"},
}


def read_cosy_column(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=1)
    return pd.to_numeric(df[0], errors="coerce").dropna()


def plot_twiss_from_dat(
    dat_dir: Path,
    *,
    title: str = "Twiss Parameters vs Element Number",
    stem: Optional[str] = None,
) -> Tuple[plt.Figure, bool]:
    fig, ax_left = plt.subplots(figsize=(12, 6))
    ax_right = ax_left.twinx()
    found = False

    for filename, cfg in FILES_TO_PLOT.items():
        full_path = dat_dir / filename
        if not full_path.is_file():
            continue
        try:
            y_values = read_cosy_column(full_path)
            ax = ax_left if cfg["ax"] == "left" else ax_right
            ax.plot(y_values.index, y_values.values, label=cfg["label"], color=cfg["color"], lw=1)
            found = True
        except Exception as exc:
            print(f"Error reading {filename}: {exc}")

    if found:
        ax_left.set_xlabel("Element number")
        ax_left.set_ylabel(r"X/Y $\beta$-function [m]")
        ax_right.set_ylabel(r"X/Y Dispersion [m]")
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        ax_left.legend(h1 + h2, l1 + l2, loc="upper right")
        ax_left.grid(True, linestyle=":", alpha=0.5)
        fig.suptitle(title if not stem else f"Twiss: {stem}")
        fig.tight_layout()

    return fig, found


def save_twiss_plot(dat_dir: Path, stem: str) -> Optional[Path]:
    fig, found = plot_twiss_from_dat(dat_dir, stem=stem)
    if not found:
        plt.close(fig)
        return None
    out_png = dat_dir / f"{stem}_twiss.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png
