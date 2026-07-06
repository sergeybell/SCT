#!/usr/bin/env python3
"""Launch COSY Infinity with cwd=COSY/src (shared runtime)."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
COSY_SRC = REPO / "COSY" / "src"
COSY_RUN = COSY_SRC / "run"
STRUCTURES = REPO / "COSY" / "structures"
DAT_ROOT = COSY_SRC / "dat"

PRE_RUN_FOX = ["cosy.fox", "utilities.fox", "elements.fox", "header.fox"]
TWISS_OUTPUT_FILES = ["BETAX", "BETAY", "DISPX"]


def _cosy_exe() -> Path:
    exe = COSY_SRC / "cosy.exe"
    if not exe.is_file():
        raise FileNotFoundError(f"cosy.exe not found at {exe}")
    return exe


def _fox_arg(path: Path) -> str:
    """Path to .fox relative to COSY/src (cosy cwd)."""
    path = path.resolve()
    try:
        return path.relative_to(COSY_SRC.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def run_cosy(fox: Path, *, verbose: bool = True) -> None:
    exe = _cosy_exe()
    arg = _fox_arg(fox)
    if verbose:
        print(f"========================================")
        print(f"RUNNING: {arg}")
        print(f"========================================")
    subprocess.run([str(exe), arg], cwd=str(COSY_SRC), check=True)


def run_pre() -> None:
    for fox_name in PRE_RUN_FOX:
        run_cosy(COSY_SRC / fox_name)


def structure_dir(stem: str) -> Path:
    return STRUCTURES / stem


def maps_fox(stem: str) -> Path:
    p = structure_dir(stem) / f"{stem}_maps.fox"
    if not p.is_file():
        raise FileNotFoundError(f"Maps file not found: {p}")
    return p


def collect_twiss_output(stem: str) -> Path:
    dest = DAT_ROOT / stem
    dest.mkdir(parents=True, exist_ok=True)
    for name in TWISS_OUTPUT_FILES:
        src = COSY_SRC / name
        if src.is_file():
            shutil.move(str(src), str(dest / name))
    return dest


def run_twiss(stem: str) -> Path:
    sys.path.insert(0, str(REPO / "COSY" / "analysis"))
    from twiss_from_fox import write_twiss_run_fox

    maps = maps_fox(stem)
    run_cosy(maps)
    twiss_fox = write_twiss_run_fox(stem)
    run_cosy(twiss_fox)
    return collect_twiss_output(stem)


def main() -> int:
    ap = argparse.ArgumentParser(description="Run COSY Infinity (cwd always COSY/src)")
    ap.add_argument("--pre", action="store_true", help="Compile base modules (cosy, utilities, elements, header)")
    ap.add_argument("--twiss", metavar="STEM", help="Run maps + Twiss for structure stem (e.g. magnetic_2)")
    ap.add_argument("fox_files", nargs="*", help=".fox paths (absolute or relative to repo)")
    args = ap.parse_args()

    if args.pre:
        run_pre()

    if args.twiss:
        dest = run_twiss(args.twiss)
        print(f"OK: Twiss output in {dest}")

    for fox in args.fox_files:
        run_cosy(Path(fox))

    if not args.pre and not args.twiss and not args.fox_files:
        ap.print_help()
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
