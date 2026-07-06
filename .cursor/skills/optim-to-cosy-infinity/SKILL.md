---
name: optim-to-cosy-infinity
description: Convert OptiM lattice files (.opt) to COSY Infinity lattice macros (.fox) in this repository. Use when the user asks to convert OptiM→COSY Infinity, mentions OptiM Magnetic lattices (magnetic_2..8.opt), or wants a quick conversion command.
disable-model-invocation: true
---

# OptiM → COSY Infinity

## Quick start

Run the converter script:

```bash
python COSY/analysis/optim_to_cosy.py OptiM/magnetic/magnetic_2.opt COSY/src/magnetic_2.fox
```

If output path is omitted, it defaults to `COSY/src/<input_stem>.fox`.

The converter writes **two** files:
- `COSY/src/<stem>.fox` — base lattice with mandatory RF block (`RFFLAG`-controlled)
- `COSY/src/<stem>_maps.fox` — same lattice with per-element `SMAPS` for Twiss (`Twiss.fox`)

## Checks

After generation, confirm the base `.fox` contains:
- `INCLUDE 'header';`
- `PROCEDURE LATTICE`
- `{SETTING RF PARAMETERS}` and `IF RFFLAG=1; RF VRF`
- `LOOP I 1 1;`
- `SAVE '<stem>';`

Confirm the maps `.fox` contains:
- `MAPARR SPNRARR` in `PROCEDURE LATTICE`
- `SMAPS 1..N` on lattice elements (RF block has no SMAPS)
- `SAVE '<stem>_maps';`
