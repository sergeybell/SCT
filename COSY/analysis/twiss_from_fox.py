#!/usr/bin/env python3
"""Build ephemeral Twiss_run.fox from TWISS SETUP block in a structure .fox file."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from optim_to_cosy import TwissSetup, _format_float

REPO = Path(__file__).resolve().parents[2]
COSY_SRC = REPO / "COSY" / "src"
TWISS_TEMPLATE = COSY_SRC / "Twiss.fox"

_RE_BLOCK_START = re.compile(r"\{---\s*TWISS SETUP")
_RE_BLOCK_END = re.compile(r"\{---\s*END TWISS\s*---\}")
_RE_KV = re.compile(r"^\{\s*([^:}]+):\s*(.+?)\s*\}$")


def _parse_kv_line(line: str) -> dict[str, str]:
    m = _RE_KV.match(line.strip())
    if not m:
        return {}
    key = m.group(1).strip()
    rest = m.group(2).strip()
    out: dict[str, str] = {}
    if ";" in rest:
        for part in rest.split(";"):
            part = part.strip()
            if not part:
                continue
            if ":" in part:
                k, v = part.split(":", 1)
                out[k.strip()] = v.strip()
            else:
                out[key] = part
    else:
        out[key] = rest
    return out


def parse_twiss_block(fox_text: str) -> TwissSetup:
    lines = fox_text.splitlines()
    start = end = None
    for i, ln in enumerate(lines):
        if _RE_BLOCK_START.search(ln):
            start = i
        if start is not None and _RE_BLOCK_END.search(ln):
            end = i
            break
    if start is None or end is None:
        raise ValueError("TWISS SETUP block not found in .fox file")

    kv: dict[str, str] = {}
    for ln in lines[start + 1 : end]:
        kv.update(_parse_kv_line(ln))

    stem = kv.get("maps_include", "").replace("_maps", "")
    if not stem:
        raise ValueError("maps_include missing in TWISS block")

    def f(key: str, default: float = 0.0) -> float:
        return float(kv.get(key, default))

    particle = kv.get("particle", "deuteron")
    set_proc = kv.get("set", "SET_FOR_DEUTERONS" if particle == "deuteron" else "SET_FOR_PROTONS")

    return TwissSetup(
        stem=stem,
        num_ele=int(kv.get("num_ele", "0")),
        maps_include=kv["maps_include"],
        particle=particle,
        gamma=f("gamma", 1.0),
        set_procedure=set_proc,
        rfflag=int(kv.get("rfflag", "0")),
        eb1=f("eb1", 0.0),
        betax0_m=f("betax0_m"),
        betay0_m=f("betay0_m"),
        alphax0=f("alphax0"),
        alphay0=f("alphay0"),
        dispx0_m=f("dispx0_m"),
        dispxpr0=f("dispxpr0"),
        dispy0_m=f("dispy0_m"),
        dispypr0=f("dispypr0"),
    )


def _extract_procedures(template_text: str) -> str:
    """Return PICK..DISPS procedures from Twiss.fox (exclude top INCLUDEs and MAIN)."""
    lines = template_text.splitlines()
    proc_start = None
    main_start = None
    for i, ln in enumerate(lines):
        if ln.strip().upper().startswith("PROCEDURE PICK"):
            proc_start = i
        if ln.strip().upper().startswith("PROCEDURE MAIN"):
            main_start = i
            break
    if proc_start is None or main_start is None:
        raise ValueError("Could not extract procedure block from Twiss.fox template")
    return "\n".join(lines[proc_start:main_start]).rstrip() + "\n"


def generate_twiss_fox(setup: TwissSetup, template_path: Optional[Path] = None) -> str:
    template = (template_path or TWISS_TEMPLATE).read_text(encoding="utf-8")
    procedures = _extract_procedures(template)

    rfflag = setup.rfflag
    main = f"""PROCEDURE MAIN;
    VARIABLE GAMMA 1;
    VARIABLE SGx1 1; VARIABLE SGx2 1; VARIABLE SGy1 1; VARIABLE SGy2 1; VARIABLE EB1 1;
    VARIABLE ALPHAX0 1; VARIABLE BETAX0 1; VARIABLE GAMMAX0 1;
    VARIABLE BETAXS 1 600;
    VARIABLE ALPHAY0 1; VARIABLE BETAY0 1; VARIABLE GAMMAY0 1;
    VARIABLE BETAYS 1 600;
    VARIABLE DISPX0 1; VARIABLE DISPXPR0 1;
    VARIABLE DISPXS 1 600;
    VARIABLE DISPY0 1; VARIABLE DISPYPR0 1;
    VARIABLE DISPYS 1 600;
    VARIABLE MAPARR1 1000 6 600; VARIABLE SPNRARR1 1000 3 3 600;
    VARIABLE NUM_ELE 1;

    GAMMA := {_format_float(setup.gamma)};
    NUM_ELE := {setup.num_ele};

    OV 3 3 0;
    {setup.set_procedure} GAMMA;

    EB1 := {_format_float(setup.eb1)};
    SGx1 := {_format_float(setup.sext_gx1)}; SGx2 := {_format_float(setup.sext_gx2)};
    SGy1 := {_format_float(setup.sext_gy1)}; SGy2 := {_format_float(setup.sext_gy2)};

    BETAX0 := {_format_float(setup.betax0_m)}; ALPHAX0 := {_format_float(setup.alphax0)};
    BETAY0 := {_format_float(setup.betay0_m)}; ALPHAY0 := {_format_float(setup.alphay0)};
    DISPX0 := {_format_float(setup.dispx0_m)}; DISPXPR0 := {_format_float(setup.dispxpr0)};
    DISPY0 := {_format_float(setup.dispy0_m)}; DISPYPR0 := {_format_float(setup.dispypr0)};

    GAMMAX0 := (1+SQR(ALPHAX0))/BETAX0;
    GAMMAY0 := (1+SQR(ALPHAY0))/BETAY0;

    LATTICE SGx1 SGy1 SGx2 SGy2 EB1 {rfflag} MAPARR1 SPNRARR1; PM 6;

    BETS 'BETAX' 'X' NUM_ELE ALPHAX0 BETAX0 GAMMAX0 BETAXS MAPARR1 SPNRARR1;
    BETS 'BETAY' 'Y' NUM_ELE ALPHAY0 BETAY0 GAMMAY0 BETAYS MAPARR1 SPNRARR1;
    DISPS 'DISPX' 'X' NUM_ELE DISPX0 DISPXPR0 DISPXS MAPARR1 SPNRARR1;

ENDPROCEDURE;

PROCEDURE RUN;
  MAIN;
ENDPROCEDURE;
RUN; END;
"""
    return (
        f"{{ Auto-generated Twiss for {setup.stem} }}\n"
        f"INCLUDE '{setup.maps_include}';\n\n"
        + procedures
        + "\n"
        + main
    )


def write_twiss_run_fox(stem: str, structures_dir: Optional[Path] = None) -> Path:
    base = structures_dir or (REPO / "COSY" / "structures" / stem)
    fox_candidates = [base / f"{stem}.fox", base / f"{stem}_maps.fox"]
    fox_path = next((p for p in fox_candidates if p.is_file()), None)
    if fox_path is None:
        raise FileNotFoundError(f"No .fox with TWISS block for stem {stem!r} in {base}")

    setup = parse_twiss_block(fox_path.read_text(encoding="utf-8"))
    out_dir = base / "_generated"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "Twiss_run.fox"
    out_path.write_text(generate_twiss_fox(setup), encoding="utf-8")
    return out_path
