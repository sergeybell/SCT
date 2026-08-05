# `py/` — пайплайн анализа

Все скрипты запускаются из корня репозитория:

```bash
python COSY/sct_study/py/<script>.py
```

Общие пути и формулы — в [`common.py`](common.py).

## Offline (без `cosy.exe`)

| Скрипт | Что делает | Основные выходы |
|--------|------------|-----------------|
| `audit_mapping.py` | Качество линейной модели mapping FR0 | `dat/<stem>/audit.json`, `plots/*_grid_vs_istar.png` |
| `validate_coefficients.py` | Конвенции η₁, MCM, γG, fACCLEN | `dat/coefficient_conventions.json` |
| `analyze_zero_point.py` | \(I^*\) и контрольные токи | `dat/<stem>/zero_point.json` |
| `analyze_delta_eq.py` | Теория \(\Delta\delta_{eq}\) и map-\(\Delta\nu_s\) | панели JSON + PNG |
| `analyze_phase_advance.py` | Полярные диаграммы секступолей + результирующий вектор | `sext_phase.json`, `*_sext_phase_{x,y}.png` |
| `analyze_working_point.py` | Парсинг `wp_*.dat`, резонансы | `working_point.json` |
| `analyze_tracking.py` | Postprocess TRPRAY/TRPSPI → \(C(n)\), \(\overline D\) и диагностический fit фазы | `track_*_analysis.json` |
| `plot_results.py` | Сводные графики | `plots/summary_*.png` |
| `run_offline_analysis.py` | Запуск всей offline-цепочки подряд | — |

## Генерация FOX и запуск COSY

| Скрипт | Что делает |
|--------|------------|
| `generate_fox.py` | Заполняет `_template_*.fox` токами и путями |
| `run_cosy_jobs.py` | `--job validate \| working_point \| track` (+ `--smoke` / `--full`) |

## Зависимости

- Python 3 + `numpy`, `pandas`, `matplotlib`
- Для COSY-jobs: `COSY/src/cosy.exe` и подготовка через `python COSY/src/run/run_cosy.py --pre`
- Mapping API: `COSY/analysis/map_lat_lib.py`

Внимание: поле `dnu_s` в существующем tracking-анализе получено по редко
сохранённой абсолютной фазе. Из-за алиасинга оно не является надёжным
индивидуальным spin tune; подробное объяснение дано в `../analysis.ipynb`.
