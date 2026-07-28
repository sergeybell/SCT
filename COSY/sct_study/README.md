# SCT study (`magnetic_2`…`magnetic_5`, FR0)

Изолированное исследование связи хроматичности, удлинения орбиты и спиновой когерентности в COSY Infinity.

## Структура

```
sct_study/
  README.md          — этот файл
  WORKLOG.md         — хронологический журнал шагов
  manifest.json      — связь stem ↔ fox/dat/plots
  config/            — параметры анализа
  fox/               — FOX для WP, I*-check, INJECT+TR
  py/                — анализ и графики
  dat/<stem>/        — результаты нового исследования
  plots/             — графики
```

Исходные mapping-данные **не копируются**: читаются из `COSY/src/dat/<stem>/`.
COSY всегда запускается с `cwd=COSY/src`.

## Быстрый старт

```bash
# 1) Audit + I* + Δδ_eq + панели (без COSY)
python COSY/sct_study/py/run_offline_analysis.py

# 2) Подготовить FOX (токи I* и контрольные точки)
python COSY/sct_study/py/generate_fox.py

# 3) COSY: валидация I* / рабочая точка / трекинг (нужен cosy.exe)
python COSY/src/run/run_cosy.py --pre
python COSY/sct_study/py/run_cosy_jobs.py --job validate
python COSY/sct_study/py/run_cosy_jobs.py --job working_point
python COSY/sct_study/py/run_cosy_jobs.py --job track --smoke

# 4) Postprocess трекинга
python COSY/sct_study/py/analyze_tracking.py
```

## Зафиксированные решения

- Только FR0 (без fringe).
- Без FIT.
- Structures: `magnetic_2`…`magnetic_5`.
- Пилот трекинга: `magnetic_2` (лучший condition number \(R_{\mathrm{int}}\)).
- Multi-turn: ваш `INJECT+TR`; анализ **индивидуальных** \(\overline D_i-\overline D_{\mathrm{ref}}\), \(\Delta\nu_{s,i}\), \(C(n)\).

См. `WORKLOG.md` для статусов и численных результатов.
