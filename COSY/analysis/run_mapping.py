#!/usr/bin/env python3
"""Generate and run mapping.fox scans per structure (FR0 / FR3) into dat/<stem>/."""
from __future__ import annotations

import argparse
import re
import shutil
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
COSY_SRC = REPO / "COSY" / "src"
STRUCTURES = REPO / "COSY" / "structures"
DAT_ROOT = COSY_SRC / "dat"
MAPPING_TEMPLATE = COSY_SRC / "mapping.fox"
MAPPING_RUN_DIR = COSY_SRC / "mapping_run"

RAW_OUTPUTS = {
    "integrals.dat": "integrals",
    "particle_spin_tune.dat": "particle_spin_tune",
    "particle_spin_tune_1.dat": "particle_spin_tune_1",
}

sys.path.insert(0, str(REPO / "COSY" / "src" / "run"))
from run_cosy import run_cosy, run_pre  # noqa: E402


def _lattice_fox(stem: str) -> Path:
    p = STRUCTURES / stem / f"{stem}.fox"
    if not p.is_file():
        raise FileNotFoundError(f"Lattice not found: {p}")
    return p


def generate_mapping_fox(stem: str, mode: str) -> str:
    if mode not in ("FR0", "FR3"):
        raise ValueError(f"mode must be FR0 or FR3, got {mode!r}")

    lines = MAPPING_TEMPLATE.read_text(encoding="utf-8").splitlines()
    proc_start = next(i for i, ln in enumerate(lines) if ln.strip().upper().startswith("PROCEDURE "))
    body = "\n".join([f"INCLUDE '{stem}';"] + lines[proc_start:]) + "\n"

    body = re.sub(
        r"^\s*structure\s*:=\s*'[^']*'\s*;",
        f"    structure := '{stem}';",
        body,
        count=1,
        flags=re.MULTILINE,
    )
    body = re.sub(r"^\s*\{structure\s*:=\s*'[^']*'\s*;\}\s*$", "", body, flags=re.MULTILINE)

    if mode == "FR0":
        body = re.sub(
            r"^\s*N0\s*:=\s*3;\s*N1\s*:=\s*3;\s*N2\s*:=\s*3;\s*\{\s*fringe field mode\s*\}",
            "{N0 :=  3; N1 := 3; N2 := 3;} { fringe field mode }",
            body,
            flags=re.MULTILINE,
        )
        fringe_line = "UM ; {FR 3 ;} LATTICE SGx1 SGy1 SGx2 SGy2 EB1 RF ; {FR 0 ;}"
    else:
        body = re.sub(
            r"N0\s*:=\s*5;\s*N1\s*:=\s*5;\s*N2\s*:=\s*5;\s*\{\s*casual mode\s*\}",
            "{N0 :=  5; N1 := 5; N2 := 5;} { casual mode }",
            body,
        )
        body = re.sub(
            r"^\s*\{N0\s*:=\s*3;\s*N1\s*:=\s*3;\s*N2\s*:=\s*3;\}\s*\{\s*fringe field mode\s*\}",
            "N0 :=  3; N1 := 3; N2 := 3; { fringe field mode }",
            body,
            flags=re.MULTILINE,
        )
        fringe_line = "UM ; FR 3 ; LATTICE SGx1 SGy1 SGx2 SGy2 EB1 RF ; FR 0 ;"

    body = body.replace(
        "UM ; {FR 3 ;} LATTICE SGx1 SGy1 SGx2 SGy2 EB1 RF ; {FR 0 ;}",
        fringe_line,
    )
    if f"structure := '{stem}'" not in body:
        raise RuntimeError(f"Failed to set structure := '{stem}'")
    return body


def _dat_dir(stem: str) -> Path:
    return DAT_ROOT / stem


def _remove_stale_raw(stem: str) -> None:
    for raw in RAW_OUTPUTS:
        p = _dat_dir(stem) / raw
        if p.is_file():
            p.unlink()


def _rename_outputs(stem: str, mode: str) -> None:
    d = _dat_dir(stem)
    d.mkdir(parents=True, exist_ok=True)
    for raw, prefix in RAW_OUTPUTS.items():
        src = d / raw
        dst = d / f"{prefix}_{mode}.dat"
        if not src.is_file():
            raise FileNotFoundError(f"Expected mapping output missing: {src}")
        if dst.is_file():
            dst.unlink()
        shutil.move(str(src), str(dst))
        print(f"  OK: {dst.relative_to(REPO)}")


def compile_lattice(stem: str) -> None:
    run_cosy(_lattice_fox(stem))


def run_mapping(stem: str, mode: str, *, skip_compile: bool = False) -> Path:
    MAPPING_RUN_DIR.mkdir(parents=True, exist_ok=True)
    fox_path = MAPPING_RUN_DIR / f"{stem}_{mode}.fox"
    fox_path.write_text(generate_mapping_fox(stem, mode), encoding="utf-8")

    _dat_dir(stem).mkdir(parents=True, exist_ok=True)
    _remove_stale_raw(stem)

    t0 = time.perf_counter()
    if not skip_compile:
        print(f"Compile lattice: {stem}")
        compile_lattice(stem)

    print(f"Mapping {stem} mode={mode}")
    run_cosy(fox_path)
    _rename_outputs(stem, mode)
    print(f"Done {stem} {mode} in {time.perf_counter() - t0:.1f}s")
    return _dat_dir(stem)


def run_mapping_modes(stem: str, modes: list[str], *, skip_compile: bool = False) -> None:
    for i, mode in enumerate(modes):
        run_mapping(stem, mode, skip_compile=skip_compile or i > 0)


def main() -> int:
    ap = argparse.ArgumentParser(description="Run COSY mapping scan for a magnetic structure")
    ap.add_argument("--stem", action="append", dest="stems", metavar="STEM")
    ap.add_argument("--mode", choices=["FR0", "FR3", "both"], default="both")
    ap.add_argument("--all-magnetic", action="store_true", help="magnetic_3..magnetic_8")
    ap.add_argument("--pre", action="store_true")
    ap.add_argument("--skip-compile", action="store_true")
    args = ap.parse_args()

    stems: list[str] = list(args.stems or [])
    if args.all_magnetic:
        stems.extend(f"magnetic_{n}" for n in range(3, 9))
    if not stems:
        ap.print_help()
        return 1

    modes = ["FR0", "FR3"] if args.mode == "both" else [args.mode]
    if args.pre:
        run_pre()

    for stem in stems:
        print(f"\n========== {stem} ==========")
        run_mapping_modes(stem, modes, skip_compile=args.skip_compile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
