#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class OptimBeam:
    energy_mev: float  # kinetic energy [MeV]
    mass_mev: float  # rest mass [MeV]

    @property
    def p_mev_c(self) -> float:
        # OptiM help: p = sqrt(E^2 + 2 E M) for kinetic energy E and rest mass M
        return math.sqrt(self.energy_mev * self.energy_mev + 2.0 * self.energy_mev * self.mass_mev)

    @property
    def brho_tm(self) -> float:
        # B*rho [T*m] = p[MeV/c] / 299.792458
        return self.p_mev_c / 299.792458


@dataclass
class OptimElement:
    name: str
    # Raw OptiM units (may be None if absent)
    l_cm: float = 0.0
    b_kg: float = 0.0
    g_kg_per_cm: float = 0.0
    s_kg_per_cm2: float = 0.0
    tilt_deg: float = 0.0

    def type_code(self) -> str:
        return (self.name[:1] or "").upper()


def _strip_comment(line: str) -> str:
    # In OptiM: lines beginning with '#' are comments.
    # We also treat trailing '#' as comment marker in math header, but we are not evaluating full header here.
    return line.rstrip("\n")


_re_assign = re.compile(r"^\s*\$(?P<name>[A-Za-z_]\w*)\s*=\s*(?P<expr>[^;#]+)")


def _parse_local_vars(lines: List[str], optim_idx: int) -> Dict[str, float]:
    """
    Parse $VAR assignments near the OptiM marker. We intentionally do NOT execute/interpret the whole math header.
    We only need simple scalars used in the element list (e.g. $QUAD).
    """
    start = max(0, optim_idx - 80)
    vars_: Dict[str, float] = {}
    for raw in lines[start : optim_idx + 40]:
        m = _re_assign.match(raw)
        if not m:
            continue
        name = m.group("name")
        expr = m.group("expr").strip()
        try:
            vars_[name] = _eval_optim_expr(expr, vars_)
        except Exception:
            # Ignore anything we can't safely evaluate; it may be notebook-like content.
            continue
    return vars_


def _eval_optim_expr(expr: str, vars_: Dict[str, float]) -> float:
    """
    Minimal evaluator for OptiM expressions used in element lists.
    Supports:
      - numbers (incl scientific notation)
      - $VAR substitution (with optional unary +/-)
      - + - * / and parentheses
      - sqrt/sin/cos/tan/abs (as in help), and pi

    This is intentionally conservative.
    """
    s = expr.strip()
    # Replace OptiM power operator '^' with Python '**'
    s = s.replace("^", "**")

    # Replace $PI/$pi and similar constants (OptiM help: $pi and $PI exist)
    s = re.sub(r"\$(PI|pi)\b", "pi", s)

    # Replace $VARS with a dictionary lookup
    def repl_var(m: re.Match) -> str:
        v = m.group(1)
        if v in ("PI", "pi"):
            return "pi"
        return f"vars_['{v}']"

    s = re.sub(r"\$([A-Za-z_]\w*)", repl_var, s)

    allowed = {
        "sqrt": math.sqrt,
        "sin": math.sin,
        "cos": math.cos,
        "tan": math.tan,
        "asin": math.asin,
        "acos": math.acos,
        "atan": math.atan,
        "exp": math.exp,
        "log": math.log,
        "abs": abs,
        "pi": math.pi,
        "vars_": vars_,
    }
    return float(eval(s, {"__builtins__": {}}, allowed))


_re_optim_marker = re.compile(r"^\s*OptiM\s*$")
_re_energy_mass = re.compile(r"Energy\[MeV\]\s*=\s*(?P<E>[-+0-9.eE$]+)\s+Mass\[MeV\]\s*=\s*(?P<M>[-+0-9.eE$]+)")


def _find_line_idx(lines: List[str], pattern: re.Pattern[str]) -> int:
    for i, ln in enumerate(lines):
        if pattern.search(ln):
            return i
    raise ValueError(f"Could not find pattern: {pattern.pattern}")


def _tokenize_lattice(lines: List[str], start_idx: int, end_idx: int) -> List[str]:
    tokens: List[str] = []
    for ln in lines[start_idx:end_idx]:
        if ln.lstrip().startswith("#"):
            continue
        # Strip any inline comment after '#'
        ln = ln.split("#", 1)[0]
        tokens.extend([t for t in re.split(r"[\s\t]+", ln.strip()) if t])
    return tokens


def _parse_element_line(line: str, vars_: Dict[str, float]) -> Optional[OptimElement]:
    if not line.strip():
        return None
    if line.lstrip().startswith("#"):
        return None
    # Example: QF1  L[cm]=21  G[kG/cm]=$QUAD  Tilt[deg]=0
    # We parse by regex key=value where value can be number or expression with $VAR.
    parts = line.split()
    if not parts:
        return None
    name = parts[0]
    rest = " ".join(parts[1:])

    def get_val(keys: Iterable[str]) -> Optional[float]:
        for k in keys:
            m = re.search(rf"{re.escape(k)}\s*=\s*([^\s]+)", rest)
            if m:
                return _eval_optim_expr(m.group(1), vars_)
        return None

    l_cm = get_val(["L[cm]"]) or 0.0
    b_kg = get_val(["B[kG]"]) or 0.0
    g_kg_per_cm = get_val(["G[kG/cm]", "Gb[kG/cm]"]) or 0.0
    s_kg_per_cm2 = get_val(["S[kG/cm/cm)]", "S[kG/cm/cm)]=", "S[kG/cm/cm)]]"]) or 0.0
    # sextupole key in files is literally "S[kG/cm/cm)]"
    if s_kg_per_cm2 == 0.0:
        mS = re.search(r"S\[kG/cm/cm\)\]\s*=\s*([^\s]+)", rest)
        if mS:
            s_kg_per_cm2 = _eval_optim_expr(mS.group(1), vars_)

    tilt_deg = get_val(["Tilt[deg]", "T[deg]"]) or 0.0

    return OptimElement(
        name=name,
        l_cm=l_cm,
        b_kg=b_kg,
        g_kg_per_cm=g_kg_per_cm,
        s_kg_per_cm2=s_kg_per_cm2,
        tilt_deg=tilt_deg,
    )


def parse_optim(path: Path) -> Tuple[OptimBeam, Dict[str, OptimElement], List[str]]:
    raw = path.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = [_strip_comment(ln) for ln in raw]

    optim_idx = _find_line_idx(lines, _re_optim_marker)
    vars_ = _parse_local_vars(lines, optim_idx)

    # Beam line is right after OptiM marker in our files, but we'll search forward within 5 lines.
    beam_line = None
    for ln in lines[optim_idx : optim_idx + 10]:
        if "Energy[MeV]" in ln and "Mass[MeV]" in ln:
            beam_line = ln
            break
    if beam_line is None:
        raise ValueError("Could not find Energy[MeV]/Mass[MeV] line after OptiM marker")

    m = _re_energy_mass.search(beam_line)
    if not m:
        raise ValueError(f"Could not parse beam line: {beam_line}")

    energy = _eval_optim_expr(m.group("E"), vars_)
    mass = _eval_optim_expr(m.group("M"), vars_)
    beam = OptimBeam(energy_mev=energy, mass_mev=mass)

    # Lattice block
    lat_begin = _find_line_idx(lines, re.compile(r"^\s*begin lattice", re.IGNORECASE))
    lat_end = _find_line_idx(lines, re.compile(r"^\s*end lattice", re.IGNORECASE))
    sequence = _tokenize_lattice(lines, lat_begin + 1, lat_end)

    # Elements list block
    list_begin = _find_line_idx(lines, re.compile(r"^\s*begin list", re.IGNORECASE))
    # end list can be "end list" or "end list of elements"
    list_end = _find_line_idx(lines, re.compile(r"^\s*end list", re.IGNORECASE))

    elements: Dict[str, OptimElement] = {}
    for ln in lines[list_begin + 1 : list_end]:
        el = _parse_element_line(ln, vars_)
        if el is None:
            continue
        elements[el.name] = el

    return beam, elements, sequence


def _cosy_header(input_path: Path, elem_count: int) -> List[str]:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return [
        "{ " + "*" * 60 + " }",
        "{ COSY Infinity Lattice Converter }",
        "{ Purpose: Conversion from OptiM to COSY Infinity }",
        "{ Author: Kolokolchikov Sergey }",
        f"{{ Date: {now} }}",
        f"{{ Input file: {input_path.as_posix()} }}",
        f"{{ TOTAL ELEMENTS: {elem_count} }}",
        "{ " + "*" * 60 + " }",
        "",
    ]


def _format_float(x: float) -> str:
    # keep compact but stable
    if abs(x) < 1e-12:
        return "0"
    return f"{x:.10g}"


_ELEMENT_KEYWORDS = ["QUAD", "SBEND", "DL", "MH", "MS", "MQ", "WIEN2D", "WIEN", "EQ", "ED", "EH"]
_RE_LATTICE_ELEMENT = re.compile(
    r"^([ \t]*)(" + "|".join(_ELEMENT_KEYWORDS) + r")\b([^;]+;)(.*)",
    re.IGNORECASE,
)
_RF_SKIP_PATTERNS = [
    re.compile(r"\{SETTING RF PARAMETERS\}", re.I),
    re.compile(r"^\s*HNUM\s*:=", re.I),
    re.compile(r"^\s*VRF\s*\(", re.I),
    re.compile(r"^\s*FREQ\s*:=", re.I),
    re.compile(r"UM;\s*CR", re.I),
    re.compile(r"IF\s+RFFLAG", re.I),
    re.compile(r"^\s*ENDIF", re.I),
]


def _is_rf_setup_line(line: str) -> bool:
    return any(p.search(line) for p in _RF_SKIP_PATTERNS)


def _rf_setup_lines() -> List[str]:
    """Mandatory RF block (pattern from magnetic_2p.fox); not included in SMAPS."""
    return [
        " {SETTING RF PARAMETERS}",
        " HNUM := 1;",
        " VRF(1, 1) := 5; {RF Voltage [kV]}",
        " FREQ := HNUM*REVFREQ(fACCLEN(1)); {RF Frequency}",
        "",
        " UM; CR;",
        " IF RFFLAG=1; RF VRF 0 FREQ 0 0.05;",
        " ENDIF;",
        "",
    ]


def generate_cosy_fox(
    beam: OptimBeam,
    elements: Dict[str, OptimElement],
    sequence: List[str],
    output_stem: str,
    input_path: Path,
) -> str:
    brho = beam.brho_tm
    A = 0.05

    # Precompute derived values for known element names.
    # For now we treat identical names as identical (OptiM tags via '&' are not present in magnetic_2.opt).
    derived: Dict[str, Dict[str, float]] = {}
    for name, el in elements.items():
        t = el.type_code()
        d: Dict[str, float] = {}
        d["L_m"] = el.l_cm / 100.0
        d["Tilt_rad"] = el.tilt_deg * math.pi / 180.0
        if t in ("B", "D"):
            B_T = el.b_kg / 10.0
            d["B_T"] = B_T
            d["ANG_rad"] = (d["L_m"] * B_T / brho) if brho != 0 else 0.0
        if t == "Q":
            # K1 [1/m^2] = (G[kG/cm] * 10 [T/m]) / Brho[T*m]
            d["K1"] = (el.g_kg_per_cm * 10.0 / brho) if brho != 0 else 0.0
        if t == "S":
            # NOTE: OptiM S is in kG/cm^2. We map to a pole-tip field estimate at radius A:
            #   B'' [T/m^2] = S[kG/cm^2] * 1000
            #   By(A) ~= 0.5 * B'' * A^2  => Bpt [T]
            # This is a best-effort physical mapping; validate on a non-zero sext lattice later.
            bpp = el.s_kg_per_cm2 * 1000.0  # T/m^2
            d["Bpt_T"] = 0.5 * bpp * (A * A)
        derived[name] = d

    out: List[str] = []
    out += _cosy_header(input_path, len(sequence))
    out += ["INCLUDE 'header';", ""]
    out += ["PROCEDURE LATTICE SEXTGx1 SEXTGy1 SEXTGx2 SEXTGy2 EB1 RFFLAG;"]

    # All VARIABLE declarations first (no assignments until block is complete).
    var_lines: List[str] = [
        " VARIABLE I 1;",
        " VARIABLE A 1;",
        " VARIABLE VRF 1 1 1;",
        " VARIABLE FREQ 1;",
        " VARIABLE HNUM 1;",
    ]
    for name in sorted(elements.keys()):
        var_lines.append(f" VARIABLE L_{name} 1;")
        t = elements[name].type_code()
        if t in ("B", "D"):
            var_lines.append(f" VARIABLE ANG_{name} 1;")
        elif t == "Q":
            var_lines.append(f" VARIABLE {name} 1;")
        elif t == "S":
            var_lines.append(f" VARIABLE {name} 1;")
    out += var_lines
    out += [""]

    out += [" {LATTICE PARAMETERS}"]
    out += [f" A := {_format_float(A)};"]
    out += [""]

    out += [" {ELEMENT VALUES FROM OPTIM}"]
    out += [
        f" {{ Beam: E={_format_float(beam.energy_mev)} MeV, M={_format_float(beam.mass_mev)} MeV, Brho={_format_float(brho)} T*m }}"
    ]
    for name in sorted(elements.keys()):
        d = derived[name]
        out += [f" L_{name} := {_format_float(d['L_m'])};"]
        t = elements[name].type_code()
        if t in ("B", "D"):
            out += [f" ANG_{name} := {_format_float(d.get('ANG_rad', 0.0))};"]
        elif t == "Q":
            out += [f" {name} := {_format_float(d.get('K1', 0.0))};"]
        elif t == "S":
            out += [f" {name} := {_format_float(d.get('Bpt_T', 0.0))};"]
    out += [""]
    out += _rf_setup_lines()

    out += [" {BEGIN LATTICE}"]
    out += [" LOOP I 1 1;"]
    for tok in sequence:
        if tok not in elements:
            out += [f"  {{! Undefined element in list: {tok} }}"]
            continue
        el = elements[tok]
        t = el.type_code()
        if t in ("O", "I"):
            out += [f"  DL L_{tok} ; {{{tok}}}"]
        elif t in ("B", "D"):
            out += [f"  SBEND L_{tok} ANG_{tok} 0 0 0 0 0 0.5 0.5 ; {{{tok}}}"]
        elif t == "Q":
            out += [f"  QUAD L_{tok} A {tok} ; {{{tok}}}"]
        elif t == "S":
            out += [f"  MH L_{tok} {tok} A ; {{{tok}}}"]
        else:
            out += [f"  {{! Unsupported element type {t} for {tok} }}"]
    out += [" ENDLOOP;"]
    out += ["ENDPROCEDURE;"]
    out += [f"SAVE '{output_stem}';"]
    out += [""]

    return "\n".join(out)


def generate_cosy_maps_fox(
    base_fox_text: str,
    *,
    input_path: Path,
    output_stem: str,
    maps_output_path: Path,
) -> str:
    """
    Transform a base .fox lattice into a maps version with per-element SMAPS.
    RF setup lines are preserved without SMAPS numbering (Twiss indexing).
    """
    lines = base_fox_text.splitlines()
    new_save_name = f"{output_stem}_maps"
    global_idx = 0
    section_counts: List[int] = []
    current_section_count = 0
    processed_lines: List[str] = []

    for line in lines:
        if "PROCEDURE LATTICE" in line.upper() and "MAPARR" not in line.upper():
            processed_lines.append(
                re.sub(
                    r"(PROCEDURE\s+LATTICE.*?)(?=\s*;)",
                    r"\1 MAPARR SPNRARR",
                    line,
                    flags=re.IGNORECASE,
                )
            )
            continue

        if "SAVE '" in line.upper() and "';" in line:
            processed_lines.append(f"SAVE '{new_save_name}';")
            continue

        section_match = re.match(r"^([ \t]*)(\{={2,}.*?={2,}.*?\})", line)
        if section_match:
            if current_section_count > 0:
                section_counts.append(current_section_count)
                current_section_count = 0
            indent = section_match.group(1)
            raw_title = section_match.group(2)
            name_match = re.search(r"={2,}([A-Za-z\s]+)", raw_title)
            name = name_match.group(1).strip() if name_match else "SECTION"
            processed_lines.append(f"{indent}{{========{name}========== elements: REPLACE_ME }}")
            continue

        if _is_rf_setup_line(line):
            processed_lines.append(line)
            continue

        elem_match = _RE_LATTICE_ELEMENT.match(line)
        if elem_match:
            global_idx += 1
            current_section_count += 1
            indent = elem_match.group(1)
            elem_cmd = elem_match.group(2) + elem_match.group(3)
            comment = elem_match.group(4).strip()
            processed_lines.append(
                f"{indent}UM; {elem_cmd} {comment} SMAPS {global_idx} MAPARR SPNRARR;"
            )
        else:
            processed_lines.append(line)

    if current_section_count > 0:
        section_counts.append(current_section_count)

    final_content = "\n".join(processed_lines) + "\n"

    final_content = final_content.replace(
        "{ Purpose: Conversion from OptiM to COSY Infinity }",
        "{ Purpose: Auto-generation of maps (SMAPS) for Twiss functions }",
    )
    final_content = re.sub(
        r"(\{ Input file: [^}]+\})\n",
        rf"\1\n{{ Output file: {maps_output_path.as_posix()} }}\n",
        final_content,
        count=1,
    )

    for count in section_counts:
        final_content = final_content.replace("REPLACE_ME", str(count), 1)

    if section_counts:
        total_str = " + ".join(map(str, section_counts))
        final_content = re.sub(
            r"\{elem-counting:.*?\}",
            f"{{elem-counting: {total_str} = {global_idx}}}",
            final_content,
        )

    return final_content


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert OptiM .opt lattice to COSY Infinity .fox")
    ap.add_argument("input", type=Path, help="Input OptiM file (.opt)")
    ap.add_argument("output", type=Path, nargs="?", help="Output COSY file (.fox). Default: COSY/src/<stem>.fox")
    ap.add_argument("--stem", type=str, default=None, help="Override SAVE stem (default: input stem)")
    ap.add_argument(
        "--maps-output",
        type=Path,
        default=None,
        help="Output maps .fox path. Default: <output_stem>_maps.fox next to base output",
    )
    args = ap.parse_args()

    beam, elements, sequence = parse_optim(args.input)

    stem = args.stem or args.input.stem
    out_path = args.output or (Path("COSY") / "src" / f"{stem}.fox")
    maps_path = args.maps_output or (out_path.parent / f"{out_path.stem}_maps.fox")

    out_text = generate_cosy_fox(beam, elements, sequence, output_stem=stem, input_path=args.input)
    maps_text = generate_cosy_maps_fox(
        out_text,
        input_path=args.input,
        output_stem=stem,
        maps_output_path=maps_path,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(out_text, encoding="utf-8")
    maps_path.write_text(maps_text, encoding="utf-8")
    print(f"OK: wrote {out_path}")
    print(f"OK: wrote {maps_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

