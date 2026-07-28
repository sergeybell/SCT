#!/usr/bin/env python3
"""Launch COSY jobs for sct_study (cwd=COSY/src)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    COSY_SRC,
    DAT_OUT,
    FOX,
    append_worklog,
    load_config,
    stems,
    update_manifest_job,
)

sys.path.insert(0, str(COSY_SRC / "run"))
from run_cosy import run_cosy, run_pre  # noqa: E402


def compile_lattice(stem: str) -> None:
    fox = COSY_SRC.parents[0] / "structures" / stem / f"{stem}.fox"
    # structures live in COSY/structures; include path from COSY/src needs SAVE already done
    # Prefer running from structures via relative path
    run_cosy(fox)


def parse_validate(path: Path) -> dict:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k = k.strip()
        v = v.strip()
        if k == "tag":
            out["tag"] = v
        else:
            out[k] = float(v.replace("D", "E"))
    return out


def parse_wp_kv(path: Path) -> dict:
    return parse_validate(path)

def job_validate(stem_filter=None) -> None:
    for stem in stems():
        if stem_filter and stem != stem_filter:
            continue
        print(f"Compile {stem}")
        compile_lattice(stem)
        zp = json.loads((DAT_OUT / stem / "zero_point.json").read_text(encoding="utf-8"))
        pred = {
            "Istar": {
                "chrom": zp["pred_chrom_at_Istar"],
                "dnu": zp["pred_dnu_s_P0_at_Istar"],
            }
        }
        for tag in ("natural", "Istar"):
            fox = FOX / f"validate_{stem}_{tag}.fox"
            print(f"Run {fox.name}")
            run_cosy(fox)
            dat = DAT_OUT / stem / f"validate_{tag}.dat"
            if not dat.is_file():
                raise FileNotFoundError(
                    f"COSY did not produce {dat}. Check FOX compile errors in {fox}"
                )
            measured = parse_validate(dat)
            # merge into zero_point
            zp.setdefault("cosy_measured", {})[tag] = measured
            if tag == "Istar":
                zp["cosy_validation"] = {
                    "status": "run",
                    "measured": measured,
                    "predicted_chrom": pred["Istar"]["chrom"],
                    "predicted_dnu_P0": pred["Istar"]["dnu"],
                    "chrom_abs_err": [
                        abs(measured["chrom_x"] - pred["Istar"]["chrom"][0]),
                        abs(measured["chrom_y"] - pred["Istar"]["chrom"][1]),
                        abs(measured["eta1"] - pred["Istar"]["chrom"][2]),
                    ],
                }
            update_manifest_job({"job": "validate", "stem": stem, "tag": tag, "dat": str(dat)})
        (DAT_OUT / stem / "zero_point.json").write_text(json.dumps(zp, indent=2), encoding="utf-8")


def job_working_point(stem_filter=None) -> None:
    for stem in stems():
        if stem_filter and stem != stem_filter:
            continue
        compile_lattice(stem)
        zp = json.loads((DAT_OUT / stem / "zero_point.json").read_text(encoding="utf-8"))
        for pt in zp["control_points"]["points"]:
            tag = pt["tag"]
            fox = FOX / f"working_point_{stem}_{tag}.fox"
            print(f"Run {fox.name}")
            run_cosy(fox)
            update_manifest_job({"job": "working_point", "stem": stem, "tag": tag})


def job_track(smoke: bool = True, stem: str | None = None) -> None:
    cfg = load_config()
    stem = stem or cfg["tracking"]["pilot_stem"]
    compile_lattice(stem)
    kind = "smoke" if smoke else "full"
    zp = json.loads((DAT_OUT / stem / "zero_point.json").read_text(encoding="utf-8"))
    for pt in zp["control_points"]["points"]:
        tag = pt["tag"]
        fox = FOX / f"track_{stem}_{tag}_{kind}.fox"
        if not fox.is_file():
            print(f"missing {fox}")
            continue
        print(f"Run {fox.name}")
        run_cosy(fox)
        update_manifest_job({"job": "track", "stem": stem, "tag": tag, "kind": kind})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", action="store_true")
    ap.add_argument("--job", choices=["validate", "working_point", "track", "all"], required=True)
    ap.add_argument("--stem", default=None)
    ap.add_argument("--smoke", action="store_true", help="for track: short NTURN")
    ap.add_argument("--full", action="store_true", help="for track: full NTURN")
    args = ap.parse_args()

    if args.pre:
        run_pre()

    if args.job in ("validate", "all"):
        job_validate(args.stem)
    if args.job in ("working_point", "all"):
        job_working_point(args.stem)
    if args.job in ("track", "all"):
        smoke = not args.full
        if args.smoke:
            smoke = True
        job_track(smoke=smoke, stem=args.stem)

    append_worklog(
        f"## COSY jobs\n\n- **Статус:** run\n- **Команда:** `python COSY/sct_study/py/run_cosy_jobs.py --job {args.job}`\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
