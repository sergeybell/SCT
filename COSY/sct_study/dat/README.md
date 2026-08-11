# `dat/` — результаты исследования

Исходный mapping читается из `COSY/src/dat/<stem>/`. Здесь — выходы SCT-study.

## Общие

| Путь | Содержание |
|------|------------|
| `coefficient_conventions.json` | γ, G, η₁, D vs Δδ_eq, MCM/fACCLEN |
| `resonance_scan_summary.json` | сводка резонансного scan |
| `generated_fox.json` | список сгенерированных FOX |

## На каждую структуру

| Файл | Содержание |
|------|------------|
| `audit.json` | cond, I*, FR0 box, **d_min / ‖ΔI‖** |
| `zero_point.json` | I*, controls, COSY-measured chrom |
| `working_point.json` | Qx, Qy, νs |
| `resonance_scan.json` | imperfection / intrinsic / synchrotron / combined |
| `delta_eq_*.json` | теоретический \(\Delta\delta_{eq}\) |
| `sext_phase.json` | геометрический proxy фазы |
| `sext_phase_twiss.json` | \(\psi=\int ds/\beta\) на экземплярах секступолей |
| `validate_*.dat`, `wp_*.dat` | сырой COSY |

## Трекинг (`magnetic_2`)

`kind ∈ {smoke, full, dense}`:

| Суффикс | Содержание |
|---------|------------|
| `track_<tag>_<kind>_PRAY.dat` | начальные условия |
| `track_<tag>_<kind>_TRPRAY.dat` | орбита (X,A,Y,B,T,**D**) |
| `track_<tag>_<kind>_TRPSPI.dat` | спин |
| `track_<tag>_<kind>_analysis.json` | `mean_D_offset`, C(n), Δν (dense), sampling/alias flags |

Dense: \(\Delta n=2\). Full/smoke: спин алиасирован; орбитальный proxy валиден.
