#!/usr/bin/env python3
"""Plot Twiss curves from COSY dat/<stem>/ output."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "COSY" / "analysis"))

from plotter_lib import save_twiss_plot  # noqa: E402

COSY_SRC = REPO / "COSY" / "src"
DAT_ROOT = COSY_SRC / "dat"


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot Twiss from COSY/src/dat/<stem>/")
    ap.add_argument("stem", help="Structure name (e.g. magnetic_2)")
    args = ap.parse_args()

    dat_dir = DAT_ROOT / args.stem
    if not dat_dir.is_dir():
        print(f"ERROR: dat folder not found: {dat_dir}")
        return 1

    out = save_twiss_plot(dat_dir, args.stem)
    if out is None:
        print(f"ERROR: no BETAX/BETAY/DISPX data in {dat_dir}")
        return 1
    print(f"OK: wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
