# `dat/magnetic_2/`

Пилотная структура (лучший `cond(R_int)≈34.5`). Здесь есть **полный** offline-анализ **и** multi-turn трекинг.

## Offline / WP / validate

- `audit.json`, `zero_point.json`, `working_point.json`
- `delta_eq_*.json`, `sext_phase.json`
- `validate_natural.dat`, `validate_Istar.dat`
- `wp_natural.dat`, `wp_Istar.dat`, `wp_ctrl_*.dat`

## Tracking (`smoke`=200 / `full`=2000 оборотов)

Для каждого tag (`natural`, `Istar`, `ctrl_xi_x`, `ctrl_xi_y`, `ctrl_eta1`):

- `track_<tag>_<kind>_{PRAY,TRPRAY,TRPSPI}.dat`
- `track_<tag>_<kind>_analysis.json` — готовый postprocess

См. родительский [`../README.md`](../README.md) про «похожие» `.dat`.
