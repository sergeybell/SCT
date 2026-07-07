#!/usr/bin/env python3
"""Quick verification for magnetic_2.opt -> magnetic_2.fox conversion."""
from __future__ import annotations

import importlib.util
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def load_converter():
    spec = importlib.util.spec_from_file_location(
        "optim_to_cosy", REPO / "COSY/analysis/optim_to_cosy.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["optim_to_cosy"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    otc = load_converter()
    opt_path = REPO / "OptiM/magnetic/magnetic_2.opt"
    fox_path = REPO / "COSY/structures/magnetic_2/magnetic_2.fox"
    maps_path = REPO / "COSY/structures/magnetic_2/magnetic_2_maps.fox"

    beam, elems, seq = otc.parse_optim(opt_path)
    brho = beam.brho_tm
    bm = elems["BM"]
    L_m = bm.l_cm / 100.0
    B_T = bm.b_kg / 10.0
    ang = L_m * B_T / brho

    print("=== OptiM parsed ===")
    print(f"sequence length: {len(seq)}")
    print(f"Energy={beam.energy_mev} MeV, Mass={beam.mass_mev} MeV")
    print(f"p={beam.p_mev_c:.6f} MeV/c, Brho={brho:.9f} T*m")
    print(f"BM: L={bm.l_cm} cm, B={bm.b_kg} kG -> L_m={L_m}, ANG={ang:.10f} rad ({math.degrees(ang):.4f} deg)")
    print(f"QF: G={elems['QF'].g_kg_per_cm} kG/cm -> K1={elems['QF'].g_kg_per_cm*10/brho:.10f}")
    print(f"QF1: G={elems['QF1'].g_kg_per_cm} kG/cm (from $QUAD=0.1935)")
    print(f"QD1: G={elems['QD1'].g_kg_per_cm} kG/cm")

    fox = fox_path.read_text(encoding="utf-8")

    def extract(name: str) -> float:
        m = re.search(rf"^\s*{re.escape(name)} := ([0-9.eE+-]+);", fox, re.MULTILINE)
        assert m, f"missing {name} in fox"
        return float(m.group(1))

    print("\n=== magnetic_2.fox values ===")
    print(f"ANG_BM={extract('ANG_BM')}")
    print(f"QF={extract('QF')}")
    print(f"QF1={extract('QF1')}")
    print(f"QD1={extract('QD1')}")
    print(f"SBEND={fox.count('SBEND')}, QUAD={fox.count('QUAD')}, MH={fox.count('MH')}, DL={fox.count('DL')}")

    assert "{SETTING RF PARAMETERS}" in fox
    assert "UM; CR;" in fox
    assert "IF RFFLAG=1; RF VRF" in fox
    assert "VARIABLE VRF 1 1 1;" in fox
    assert "SF1 := SEXTGx1" in fox and "SF2 := SEXTGx2" in fox and "SD := SEXTGy1" in fox

    # cross-check OptiM header notebook value (not used by converter)
    notebook_br = 3.47648969  # from $BR in magnetic_2.opt line 159
    print(f"\nNotebook $BR in .opt header: {notebook_br} T*m (converter Brho: {brho:.9f})")

    smoke_path = REPO / "COSY/src/magnetic_2_smoke.fox"
    if smoke_path.exists():
        print("\nWARNING: magnetic_2_smoke.fox exists (leftover smoke-test artifact; delete it).")

    maps = maps_path.read_text(encoding="utf-8")
    smaps_indices = [int(m) for m in re.findall(r"SMAPS (\d+) MAPARR", maps)]
    print("\n=== magnetic_2_maps.fox ===")
    print(f"SMAPS count={len(smaps_indices)}, last index={smaps_indices[-1] if smaps_indices else 0}")
    assert "MAPARR SPNRARR" in maps
    assert "SAVE 'magnetic_2_maps';" in maps
    assert "IF RFFLAG=1; RF VRF" in maps
    assert "SMAPS" not in re.split(r"IF RFFLAG=1", maps, maxsplit=1)[0].split("UM; CR;")[-1]
    assert len(smaps_indices) == len(seq)
    assert smaps_indices == list(range(1, len(seq) + 1))

    # tolerances
    assert abs(extract("ANG_BM") - ang) < 1e-8
    assert abs(extract("QF") - elems["QF"].g_kg_per_cm * 10 / brho) < 1e-8
    assert len(seq) == 84
    assert fox.count("MH") == seq.count("SF1") + seq.count("SF2") + seq.count("SD")
    print("\nAll checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
