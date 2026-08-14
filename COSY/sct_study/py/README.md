# `py/` — пайплайн SCT study

| Скрипт | Роль |
|--------|------|
| `common.py` | пути, \(\Delta\delta_{eq}\), `currents_in_box` + нормированный \(d_{\min}\) |
| `audit_mapping.py` | FR0 audit + grid vs I* |
| `validate_coefficients.py` | конвенции MCM/η₁/D |
| `analyze_zero_point.py` | I* и control points |
| `analyze_delta_eq.py` | теоретические панели \(\Delta\delta_{eq}\) |
| `analyze_phase_advance.py` | геометрический proxy \(2\pi Qs/C\) |
| `analyze_twiss_phase.py` | \(\psi=\int ds/\beta\) на каждом экземпляре секступоля |
| `analyze_optim_sext_phase.py` | OptiM \(\mathrm{Nu}\) vs Twiss vs \(2\pi Q_{\mathrm{full}} s/C\); \(A\): SF1,SF2\(>0\), SD\(<0\) |
| `analyze_working_point.py` | Qx,Qy,νs из COSY WP |
| `analyze_resonances.py` | imperfection / intrinsic / synchrotron / combined |
| `generate_fox.py` | validate / WP / track (smoke, full, dense) |
| `run_cosy_jobs.py` | запуск COSY (`--dense` для спина) |
| `analyze_tracking.py` | `mean_D_offset`, C(n), \(\Delta\nu_s\) (dense) |
| `run_offline_analysis.py` | все offline шаги подряд |
| `plot_results.py` | сводные картинки |

## Важно

- `mean_D_offset` ≠ теоретический \(\Delta\delta_{eq}\).
- Spin tune только при `save_every ≤ 3` и `psi_deg≈0`.
- Twiss-фаза — основной результат; geometric polar — preview.
