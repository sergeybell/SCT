# `config/` — параметры исследования

| Файл | Назначение |
|------|------------|
| [`study_config.json`](study_config.json) | stems, γ, G, длины, трекинг, резонансы |

## Ключевые поля

- **`gamma`, `G`** — дейтрон; \(\gamma G\) задаёт спиновый tune.
- **`lattice_lengths_m`** — геометрическая \(L\) для \(\Delta\delta_{eq}\).
- **`fACCLEN_header=141`** — legacy RF length; ≠ геометрической \(L\).
- **`tracking`**:
  - `pilot_stem=magnetic_2`
  - `psi_deg=0` — горизонтальный/продольный начальный спин (нужен для `atan2(Sx,Sz)`)
  - `smoke_nturn=200`, `nturn=2000`, `dense_nturn=2000`
  - `save_every_turns_smoke=10`, `save_every_turns_sparse=100`, `save_every_turns_dense=2`
  - Найквист: \(\Delta n\le 1/(2|\nu_s|)\approx 3\)
- **`Qs_estimate.default`** — временная оценка для synchrotron sidebands.

Читается через `py/common.py` → `load_config()`.
