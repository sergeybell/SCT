# `config/` — параметры исследования

| Файл | Назначение |
|------|------------|
| [`study_config.json`](study_config.json) | Единый конфиг: stems, γ, G, длины решёток, сценарии эмиттанса, параметры трекинга |

## Ключевые поля `study_config.json`

- **`stems`** — `magnetic_2`…`magnetic_5`.
- **`gamma`, `G`** — энергия дейтрона и аномальный магнитный момент; \(\gamma G\) задаёт спиновый tune на рабочей точке.
- **`lattice_lengths_m`** — геометрическая длина кольца \(L\) (используется в \(\Delta\delta_{eq}\)).
- **`fACCLEN_header`** — значение `141` из заголовка FOX; **не равно** геометрическому \(L\) (известное предупреждение).
- **`emittance_scenarios`** — диагностические \(\varepsilon_x,\varepsilon_y,\delta_m\) для теоретических панелей \(\Delta\delta_{eq}\).
- **`tracking`**:
  - `pilot_stem`: сейчас только `magnetic_2`;
  - `num_per_group=9` → ~29 лучей (ref + X + Y + D);
  - `psi_deg=90` → начальный спин \(S=(0,1,0)\);
  - `smoke_nturn=200`, `nturn=2000`;
  - `save_every_frac=0.05` — доля оборотов между снимками.

Конфиг читается всеми скриптами через `py/common.py` → `load_config()`.
