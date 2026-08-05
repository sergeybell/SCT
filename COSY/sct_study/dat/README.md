# `dat/` — результаты исследования

Исходный mapping **не лежит здесь**: он читается из `COSY/src/dat/<stem>/`.  
Здесь — только выходы SCT-study.

## Общие файлы

| Путь | Содержание |
|------|------------|
| `coefficient_conventions.json` | γ, G, η₁, предупреждения MCM/fACCLEN |

## На каждую структуру `magnetic_2`…`magnetic_5`

| Файл | Содержание |
|------|------------|
| `audit.json` | `cond(R_int)`, линейность, I* vs grid |
| `zero_point.json` | \(I^*\), control points, (после COSY) measured chrom |
| `working_point.json` | \(Q_x,Q_y,\nu_s\), расстояния до резонансов |
| `delta_eq_*.json` | Теоретические панели \(\Delta\delta_{eq}\) / map-\(\Delta\nu_s\) |
| `sext_phase.json` | Позиции секступолей, \(\mu\), \(K_2\), результирующий вектор |
| `validate_{natural,Istar}.dat` | Сырой вывод COSY validate |
| `wp_{natural,Istar,ctrl_*}.dat` | Сырой вывод COSY working point |

## Только `magnetic_2` — multi-turn трекинг

Для каждого `tag` ∈ `{natural, Istar, ctrl_xi_x, ctrl_xi_y, ctrl_eta1}` и `kind` ∈ `{smoke, full}`:

| Суффикс | Содержание |
|---------|------------|
| `track_<tag>_<kind>_PRAY.dat` | Начальные условия ~29 лучей |
| `track_<tag>_<kind>_TRPRAY.dat` | Орбита по сохранённым оборотам |
| `track_<tag>_<kind>_TRPSPI.dat` | Спин по тем же оборотам |
| `track_<tag>_<kind>_analysis.json` | Postprocess: группы, \(\overline D\), \(\Delta\nu_{s,i}\), \(C(n)\) |

### Почему файлы кажутся одинаковыми

1. **`PRAY` smoke ≡ `PRAY` full** — стартовые координаты не зависят от `NTURN`.
2. **`TRPRAY`/`TRPSPI` smoke и full имеют почти одинаковое число строк**, потому что шаг записи масштабируется: `NINT(NTURN·0.05)` → ~20 снимков на луч и для 200, и для 2000 оборотов. **Содержимое разное**: последние `turn` равны 200 vs 2000.
3. Между `tag` одинакова **топология** ансамбля, меняются токи секступолей и динамика.
4. Между stems одинаков **формат** `wp_*` / `validate_*`, меняются числа.

## README по stems

См. также краткие README внутри `magnetic_2/`…`magnetic_5/`.
