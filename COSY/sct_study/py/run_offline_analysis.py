#!/usr/bin/env python3
"""Run all offline (no cosy.exe) analysis steps."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PY = Path(__file__).resolve().parent
SCRIPTS = [
    "audit_mapping.py",
    "validate_coefficients.py",
    "analyze_zero_point.py",
    "analyze_delta_eq.py",
    "analyze_phase_advance.py",
    "analyze_twiss_phase.py",
    "analyze_optim_sext_phase.py",
    "generate_fox.py",
    "analyze_working_point.py",
    "analyze_resonances.py",
    "plot_results.py",
    "analyze_tracking.py",
]


def main() -> int:
    for name in SCRIPTS:
        print(f"\n===== {name} =====")
        r = subprocess.run([sys.executable, str(PY / name)], check=False)
        if r.returncode != 0:
            print(f"FAILED {name} code={r.returncode}")
            return r.returncode
    print("\nAll offline analysis scripts finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
