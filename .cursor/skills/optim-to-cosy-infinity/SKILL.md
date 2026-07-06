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

## Checks

After generation, confirm the output `.fox` contains:
- `INCLUDE 'header';`
- `PROCEDURE LATTICE`
- `LOOP I 1 1;`
- `SAVE '<stem>';`

