# SCT study (`magnetic_2`…`magnetic_5`, FR0)

Изолированное исследование связи **орбитальной хроматичности**, **удлинения пути** и **спиновой когерентности** в COSY Infinity для дейтрона.

## Физическая цель

1. По mapping FR0 найти токи секступолей \(I^\*=(SGx1,SGx2,SGy1)\), при которых \(\xi_x=\xi_y=\eta_1=0\).
2. Через теорию \(\Delta\delta_{eq}\) (Senichev) оценить подавление удлинения орбиты.
3. Multi-turn трекингом (пилот `magnetic_2`) проверить proxy `mean_D_offset` и спиновые метрики на **плотной** выборке.

## Словарь (важно)

| Термин | Смысл |
|--------|--------|
| **\(D\)** | 6-я координата COSY в `TRPRAY` (относительное \(\delta p/p\)), не \(\Delta\delta_{eq}\). |
| **\(\Delta\delta_{eq}\)** | Теоретический сдвиг равновесного импульса из \(\xi,\eta_1,\varepsilon\). |
| **`mean_D_offset`** | \(\overline D_i-\overline D_{\rm ref}\) — трекинговый proxy. |
| **\(C(n)\)** | Фазовая когерентность в сохранённые моменты. |
| **\(\Delta\nu_s\)** | Относительный spin tune; только при \(\Delta n_{\rm save}\le 1/(2|\nu_s|)\approx 3\). |
| **smoke / full / dense** | 200 / 2000 оборотов; full — редкое сохранение (орбита); dense — \(\Delta n=2\) (спин). |
| **psi_deg=0** | Начальный спин в горизонтальной плоскости; `atan2(Sx,Sz)` измерим. При 90° (вертикаль) фаза бессмысленна. |

Подробный разбор: [`analysis.ipynb`](analysis.ipynb).

## Структура

```
sct_study/
  README.md, WORKLOG.md, analysis.ipynb, manifest.json
  config/   — study_config.json
  fox/      — шаблоны и сгенерированные .fox
  py/       — пайплайн
  dat/      — JSON и .dat
  plots/    — графики
```

Mapping читается из `COSY/src/dat/<stem>/`. COSY запускается с `cwd=COSY/src`.

## Быстрый старт

```bash
# Offline
python COSY/sct_study/py/run_offline_analysis.py

# FOX + COSY
python COSY/sct_study/py/generate_fox.py
python COSY/src/run/run_cosy.py --pre
python COSY/sct_study/py/run_cosy_jobs.py --job validate
python COSY/sct_study/py/run_cosy_jobs.py --job working_point
python COSY/sct_study/py/run_cosy_jobs.py --job track --dense --stem magnetic_2
python COSY/sct_study/py/analyze_tracking.py
python COSY/sct_study/py/analyze_resonances.py
```

## Зафиксированные решения

- Только FR0, без FIT; stems `magnetic_2`…`magnetic_5`; пилот `magnetic_2`.
- Spin tune только из `dense` (\(\Delta n=2\)); sparse `full` — для орбиты.
- MCM2: конвенция `mapping.fox` (без факториала 2).
- `fACCLEN=141` ≠ геометрическая \(L\) — известное предупреждение для RF/\(Q_s\).

См. [`WORKLOG.md`](WORKLOG.md).
