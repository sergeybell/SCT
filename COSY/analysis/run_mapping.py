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

# Per-structure, per-mode scan steps (SGx1/SGx2 step_x, SGy1 step_y).
# FR0: 5^3 grid (corner index 4); FR3: 3^3 grid (corner index 2) — larger steps.
# Three families only: SF1<-SEXTGx1, SF2<-SEXTGx2, SD<-SEXTGy1.
DEFAULT_STEP_X = 0.02
DEFAULT_STEP_Y = -0.03
CHROM_ABS_LIMIT = 20.0

MappingSteps = dict[str, dict[str, tuple[float, float]]]

MAPPING_STEPS: MappingSteps = {
    "magnetic_2": {
        "FR0": (0.0044, -0.0066),
        "FR3": (0.0078, -0.0117),
    },
    "magnetic_3": {
        "FR0": (0.0090, -0.0135),
        "FR3": (0.0153, -0.0230),
    },
    "magnetic_4": {
        "FR0": (0.0155, -0.0232),
        "FR3": (0.0255, -0.0383),
    },
    "magnetic_5": {
        "FR0": (0.0222, -0.0333),
        "FR3": (0.0377, -0.0566),
    },
    "magnetic_6": {
        "FR0": (0.0337, -0.0505),
        "FR3": (0.0533, -0.0800),
    },
    "magnetic_7": {
        "FR0": (0.0444, -0.0667),
        "FR3": (0.0703, -0.1055),
    },
    "magnetic_8": {
        "FR0": (0.0595, -0.0893),
        "FR3": (0.0909, -0.1363),
    },
    # Legacy COSY/src lattices (electrostatic ED, Nuclotron Wien filter)
    "electrostatic": {
        "FR0": (1.0, -3.0),
        "FR3": (1.0, -3.0),
    },
    "Nuclotron_8": {
        "FR0": (0.02, -0.03),
        "FR3": (0.04, -0.06),
    },
    "Nuclotron_16": {
        "FR0": (0.001, -0.001),
        "FR3": (0.002, -0.003),
    },
}

# EB1 passed to LATTICE (ED field / Wien E). Magnetic stems keep 0.
EB1_BY_STEM: dict[str, float] = {
    "electrostatic": 112.464392,
    "Nuclotron_8": 132.0,
    "Nuclotron_16": 132.0,
}

PRESENTATION_STEMS = [
    "magnetic_2",
    "magnetic_3",
    "electrostatic",
    "Nuclotron_8",
    "Nuclotron_16",
]

# COSY SAVE name may differ from stem (used for INCLUDE after compile).
INCLUDE_NAME: dict[str, str] = {
    "Nuclotron_16": "seq_Nuclotron_16",
}

# Particle / gamma overrides (default: deuteron, GAMMA from mapping.fox).
PARTICLE_SETUP: dict[str, dict[str, str | float]] = {
    "electrostatic": {
        "gamma": 1.24810736,  # proton ~233 MeV (Plotter / Twiss.fox)
        "chrom_set": "SET_FOR_PROTONS_CHROM",
        "spin_set": "SET_FOR_PROTONS",
    },
}

sys.path.insert(0, str(REPO / "COSY" / "src" / "run"))
from run_cosy import run_cosy, run_pre  # noqa: E402


def _lattice_fox(stem: str) -> Path:
    """Resolve lattice fox: structures/<stem>/ first, then COSY/src/<stem>.fox."""
    candidates = [
        STRUCTURES / stem / f"{stem}.fox",
        COSY_SRC / f"{stem}.fox",
    ]
    for p in candidates:
        if p.is_file():
            return p
    checked = ", ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"Lattice not found for {stem}: {checked}")


def get_mapping_steps(stem: str, mode: str) -> tuple[float, float]:
    """Return (step_x, step_y) for a structure and fringe mode."""
    if stem in MAPPING_STEPS and mode in MAPPING_STEPS[stem]:
        return MAPPING_STEPS[stem][mode]
    return DEFAULT_STEP_X, DEFAULT_STEP_Y


def get_eb1(stem: str) -> float:
    return float(EB1_BY_STEM.get(stem, 0.0))


def generate_mapping_fox(stem: str, mode: str) -> str:
    if mode not in ("FR0", "FR3"):
        raise ValueError(f"mode must be FR0 or FR3, got {mode!r}")

    lines = MAPPING_TEMPLATE.read_text(encoding="utf-8").splitlines()
    proc_start = next(i for i, ln in enumerate(lines) if ln.strip().upper().startswith("PROCEDURE "))
    include_name = INCLUDE_NAME.get(stem, stem)
    body = "\n".join([f"INCLUDE '{include_name}';"] + lines[proc_start:]) + "\n"

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

    step_x, step_y = get_mapping_steps(stem, mode)
    body = re.sub(
        r"^\s*step_x\s*:=\s*[^;]+;\s*step_y\s*:=\s*[^;]+;",
        f"    step_x := {step_x}; step_y := {step_y} ;",
        body,
        count=1,
        flags=re.MULTILINE,
    )

    eb1 = get_eb1(stem)
    body = re.sub(
        r"^\s*EB1\s*:=\s*[^;]+;",
        f"    EB1 := {eb1} ;",
        body,
        count=1,
        flags=re.MULTILINE,
    )

    setup = PARTICLE_SETUP.get(stem)
    if setup:
        body = re.sub(
            r"^\s*GAMMA\s*:=\s*[^;]+;",
            f"    GAMMA := {setup['gamma']};",
            body,
            count=1,
            flags=re.MULTILINE,
        )
        body = body.replace(
            "SET_FOR_DEUTERONS_CHROM GAMMA",
            f"{setup['chrom_set']} GAMMA",
        )
        body = body.replace(
            "SET_FOR_DEUTERONS GAMMA",
            f"{setup['spin_set']} GAMMA",
        )

    if f"structure := '{stem}'" not in body:
        raise RuntimeError(f"Failed to set structure := '{stem}'")
    if f"INCLUDE '{include_name}';" not in body.split("PROCEDURE", 1)[0] and not body.startswith(
        f"INCLUDE '{include_name}';"
    ):
        # INCLUDE is always the first line
        if not body.lstrip().startswith(f"INCLUDE '{include_name}';"):
            raise RuntimeError(f"Failed to INCLUDE '{include_name}' for {stem}")
    return body


def _validate_outputs(stem: str, mode: str, *, min_bytes: int = 500) -> None:
    """Fail fast if COSY wrote only headers (empty scan)."""
    d = _dat_dir(stem)
    for raw, prefix in RAW_OUTPUTS.items():
        # after rename the mode-suffixed file must exist and be non-trivial
        p = d / f"{prefix}_{mode}.dat"
        if not p.is_file() or p.stat().st_size < min_bytes:
            raise RuntimeError(
                f"Mapping output too small or missing for {stem} {mode}: {p} "
                f"(size={p.stat().st_size if p.is_file() else 'n/a'}). "
                "Check particle/EB1/INCLUDE settings and COSY log."
            )


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

    step_x, step_y = get_mapping_steps(stem, mode)
    print(f"Mapping {stem} mode={mode} (3-family, step_x={step_x}, step_y={step_y})")
    run_cosy(fox_path)
    _rename_outputs(stem, mode)
    _validate_outputs(stem, mode)
    print(f"Done {stem} {mode} in {time.perf_counter() - t0:.1f}s")
    return _dat_dir(stem)


def run_mapping_modes(stem: str, modes: list[str], *, skip_compile: bool = False) -> None:
    for i, mode in enumerate(modes):
        run_mapping(stem, mode, skip_compile=skip_compile or i > 0)


def main() -> int:
    ap = argparse.ArgumentParser(description="Run COSY mapping scan for a lattice structure")
    ap.add_argument("--stem", action="append", dest="stems", metavar="STEM")
    ap.add_argument("--mode", choices=["FR0", "FR3", "both"], default="both")
    ap.add_argument(
        "--all-magnetic",
        action="store_true",
        help="all magnetic_* stems that have an entry in MAPPING_STEPS",
    )
    ap.add_argument(
        "--presentation",
        action="store_true",
        help="Mapping_presentation set: magnetic_2/3, electrostatic, Nuclotron_8/16",
    )
    ap.add_argument("--pre", action="store_true")
    ap.add_argument("--skip-compile", action="store_true")
    args = ap.parse_args()

    stems: list[str] = list(args.stems or [])
    if args.presentation:
        stems.extend(PRESENTATION_STEMS)
    if args.all_magnetic:
        stems.extend(sorted(s for s in MAPPING_STEPS if s.startswith("magnetic_")))
    # de-dupe preserving order
    seen: set[str] = set()
    stems = [s for s in stems if not (s in seen or seen.add(s))]
    if not stems:
        ap.print_help()
        return 1

    modes = ["FR0", "FR3"] if args.mode == "both" else [args.mode]
    if args.pre:
        run_pre()

    for stem in stems:
        print(f"\n========== {stem} ==========")
        print(f"  lattice={_lattice_fox(stem)}")
        print(f"  EB1={get_eb1(stem)}")
        run_mapping_modes(stem, modes, skip_compile=args.skip_compile)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
