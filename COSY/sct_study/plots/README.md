# `plots/` — графики SCT study

| Шаблон имени | Содержание |
|--------------|------------|
| `<stem>_grid_vs_istar.png` | Сетка mapping vs точка \(I^*\) (3 проекции токов) |
| `<stem>_delta_eq_panels.png` | Теория \(\Delta\delta_{eq}\) vs \(\xi_x,\xi_y,\eta_1\) |
| `<stem>_map_dnu_panels.png` | Map-уровень \(\Delta\nu_s\) vs хроматичность |
| `<stem>_sext_phase_x.png` / `_y.png` | Полярные векторы секступолей (отдельный масштаб на stem) + результирующий вектор |
| `working_points_summary.png` | Сводка рабочих точек |
| `summary_delta_eq_suppression.png` | Сравнение подавления \(\Delta\delta_{eq}\) |
| `magnetic_2_track_<tag>_<kind>.png` | 4 панели: \(\Delta\delta_{eq}\), \(\Delta\nu_{s,i}\), \(C(n)\), RMS фазы |
| `magnetic_2_track_compare_<kind>.png` | Сравнение \(C(n)\) и RMS \(\Delta\delta_{eq}\) по tag |

Графики строятся скриптами из `py/` и дублируются/поясняются в корневом [`../analysis.ipynb`](../analysis.ipynb).
