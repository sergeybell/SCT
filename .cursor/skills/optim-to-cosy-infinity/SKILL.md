---
name: optim-to-cosy-infinity
description: Convert OptiM lattice files (.opt) to COSY Infinity lattice macros (.fox) in this repository. Use when the user asks to convert OptiM→COSY Infinity, mentions OptiM Magnetic lattices (magnetic_2..8.opt), or wants a quick conversion command.
disable-model-invocation: true
---

# OptiM → COSY Infinity

## Quick start

Run the converter script:

```bash
python COSY/analysis/optim_to_cosy.py OptiM/magnetic/magnetic_2.opt
```

Default output: `COSY/structures/<stem>/<stem>.fox` and `<stem>_maps.fox`.

Optional explicit output path:

```bash
python COSY/analysis/optim_to_cosy.py OptiM/magnetic/magnetic_2.opt COSY/structures/magnetic_2/magnetic_2.fox
```

The converter writes **two** files and embeds a `{--- TWISS SETUP (from OptiM) ---}` block with initial β, α, dispersion (cm→m), γ, particle type.

## Run COSY Twiss

`cosy.exe` must run with `cwd=COSY/src`. Use the launcher (does not replace `run.ipynb`):

```bash
python COSY/src/run/run_cosy.py --pre
python COSY/src/run/run_cosy.py --twiss magnetic_2
```

Twiss data lands in `COSY/src/dat/<stem>/` (BETAX, BETAY, DISPX).

Plot:

```bash
python COSY/analysis/plot_twiss.py magnetic_2
```

Saves `COSY/src/dat/<stem>/<stem>_twiss.png`.

## Layout

```
COSY/structures/<stem>/     # .fox from converter + _generated/Twiss_run.fox
COSY/src/dat/<stem>/        # Twiss output + PNG
COSY/src/                   # cosy.exe, header, legacy fox
```

## Checks

After generation, confirm the base `.fox` contains:
- `INCLUDE 'header';`
- `{--- TWISS SETUP (from OptiM) ---}`
- `PROCEDURE LATTICE`
- `{SETTING RF PARAMETERS}` and `IF RFFLAG=1; RF VRF`
- `SAVE '<stem>';`

Confirm the maps `.fox` contains:
- `MAPARR SPNRARR` in `PROCEDURE LATTICE`
- `SMAPS 1..N` on lattice elements (RF block has no SMAPS)
- `SAVE '<stem>_maps';`
