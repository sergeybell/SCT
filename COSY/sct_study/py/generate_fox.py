#!/usr/bin/env python3
"""Fill FOX templates for validate / working_point / track jobs."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DAT_OUT, FOX, append_worklog, load_config, stems, write_json  # noqa: E402


def render(template: str, mapping: dict) -> str:
    out = template
    for k, v in mapping.items():
        out = out.replace("{{" + k + "}}", str(v))
    return out


def fox_rel(path: Path) -> str:
    """Path relative to COSY/src for OPENF."""
    # COSY cwd = COSY/src → ../sct_study/...
    return "../sct_study/" + path.relative_to(path.parents[2]).as_posix()


def main() -> int:
    cfg = load_config()
    gamma = cfg["gamma"]
    tr = cfg["tracking"]
    generated = []

    tpl_val = (FOX / "_template_validate.fox").read_text(encoding="utf-8")
    tpl_wp = (FOX / "_template_working_point.fox").read_text(encoding="utf-8")
    tpl_tr = (FOX / "_template_track.fox").read_text(encoding="utf-8")

    for stem in stems():
        zp = json.loads((DAT_OUT / stem / "zero_point.json").read_text(encoding="utf-8"))
        for pt in zp["control_points"]["points"]:
            tag = pt["tag"]
            I = pt["I"]
            # validate only for Istar + natural (all stems)
            if tag in ("Istar", "natural"):
                outdat = DAT_OUT / stem / f"validate_{tag}.dat"
                text = render(
                    tpl_val,
                    {
                        "STEM": stem,
                        "GAMMA": gamma,
                        "SGX1": I[0],
                        "SGX2": I[1],
                        "SGY1": I[2],
                        "TAG": tag,
                        "OUTDAT": f"../sct_study/dat/{stem}/validate_{tag}.dat",
                    },
                )
                path = FOX / f"validate_{stem}_{tag}.fox"
                path.write_text(text, encoding="utf-8")
                generated.append(str(path))

            # working point for all control tags
            text = render(
                tpl_wp,
                {
                    "STEM": stem,
                    "GAMMA": gamma,
                    "SGX1": I[0],
                    "SGX2": I[1],
                    "SGY1": I[2],
                    "TAG": tag,
                    "OUTDAT": f"../sct_study/dat/{stem}/wp_{tag}.dat",
                },
            )
            path = FOX / f"working_point_{stem}_{tag}.fox"
            path.write_text(text, encoding="utf-8")
            generated.append(str(path))

        # tracking only for pilot stem
        if stem != tr["pilot_stem"]:
            continue
        for pt in zp["control_points"]["points"]:
            tag = pt["tag"]
            I = pt["I"]
            for nturn, kind in [(tr["smoke_nturn"], "smoke"), (tr["nturn"], "full")]:
                text = render(
                    tpl_tr,
                    {
                        "STEM": stem,
                        "GAMMA": gamma,
                        "SGX1": I[0],
                        "SGX2": I[1],
                        "SGY1": I[2],
                        "NTURN": nturn,
                        "PNUM": tr["num_per_group"],
                        "PSI": tr["psi_deg"],
                        "XAMP": tr["x_amp"],
                        "DSCALE": tr["d_scale"],
                        "SAVEFRAC": tr["save_every_frac"],
                        "OUTPRAY": f"../sct_study/dat/{stem}/track_{tag}_{kind}_PRAY.dat",
                        "OUTTRPRAY": f"../sct_study/dat/{stem}/track_{tag}_{kind}_TRPRAY.dat",
                        "OUTTRPSPI": f"../sct_study/dat/{stem}/track_{tag}_{kind}_TRPSPI.dat",
                    },
                )
                path = FOX / f"track_{stem}_{tag}_{kind}.fox"
                path.write_text(text, encoding="utf-8")
                generated.append(str(path))

    write_json(DAT_OUT / "generated_fox.json", {"files": generated})
    append_worklog(
        "## Generate FOX\n\n"
        f"- **Статус:** prepared\n"
        f"- **Команда:** `python COSY/sct_study/py/generate_fox.py`\n"
        f"- **Создано файлов:** {len(generated)}\n"
        f"- Пилот трекинга: `{tr['pilot_stem']}`\n"
    )
    print(f"OK: generated {len(generated)} fox files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
