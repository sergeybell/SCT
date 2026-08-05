# SCT study (`magnetic_2`…`magnetic_5`, FR0)

Изолированное исследование связи **орбитальной хроматичности**, **удлинения пути** и **спиновой когерентности** в COSY Infinity для дейтрона.

## Физическая цель

1. По mapping FR0 найти токи секступолей \(I^\*=(SGx1,SGx2,SGy1)\), при которых \(\xi_x=\xi_y=\eta_1=0\).
2. Через теорию \(\Delta\delta_{eq}\) оценить, как обнуление хроматичности подавляет удлинение орбиты.
3. Multi-turn трекингом (пилот `magnetic_2`) проверить proxy \(\overline{D}_i-\overline{D}_{ref}\) и спиновую метрику \(C(n)\).

**Важно:** \(C(n)\) — это **фазовая** когерентность (одинаковость мгновенной фазы спина относительно reference). Одинаковость физических spin tune — отдельная задача. Существующее поле `dnu_s` получено из редко сохранённой фазы и подвержено временному алиасингу, поэтому его нельзя считать надёжной индивидуальной частотой. Орбитальные \(Q_{x,i},Q_{y,i}\) из трекинга текущий код также не извлекает.

Подробный разбор модели и математики: [`analysis.ipynb`](analysis.ipynb).

## Структура и документация

```
sct_study/
  README.md              — этот файл
  WORKLOG.md             — журнал шагов и численные результаты
  analysis.ipynb         — единое полотно: физика + код + графики
  manifest.json          — реестр stem ↔ fox/dat/plots и журнал jobs
  config/README.md       — параметры (γ, G, smoke/full, ансамбль)
  fox/README.md          — шаблоны и сгенерированные .fox
  py/README.md           — скрипты пайплайна
  dat/README.md          — JSON и .dat результаты (+ почему файлы «похожи»)
  plots/README.md        — графики
```

Исходные mapping-данные **не копируются**: читаются из `COSY/src/dat/<stem>/`.  
COSY всегда запускается с `cwd=COSY/src`.

## Словарь

| Термин | Смысл |
|--------|--------|
| **smoke** | Короткий трекинг: `NTURN=200`. Проверка пайплайна («дымовой тест»). |
| **full** | Полный трекинг: `NTURN=2000`. |
| **natural** | Секступоли выключены: \(I=0\). |
| **Istar** | Токи \(I^\*\), обнуляющие \(\xi_x,\xi_y,\eta_1\) в линейной модели. |
| **ctrl_xi_x / ctrl_xi_y / ctrl_eta1** | Точки сетки с большим \(|\xi_x|\), \(|\xi_y|\) или \(|\eta_1|\) для сравнения. |
| **PRAY** | Начальные условия лучей (одинаковы для smoke и full). |
| **TRPRAY** | Орбита \(X,A,Y,B,T,D\) по сохранённым оборотам. |
| **TRPSPI** | Спин \(S_x,S_y,S_z\) по тем же оборотам. |
| **C(n)** | \(\bigl\|N^{-1}\sum_i e^{i(\phi_i-\phi_{ref})}\bigr\|\in[0,1]\). |

Почему smoke/full `.dat` выглядят «одинаковыми»: шаг сохранения `NINT(NTURN·0.05)` → ~одинаковое число строк, но разные номера оборотов (до 200 vs 2000). `PRAY` действительно идентичен.

## Быстрый старт

```bash
# 1) Offline (без cosy.exe): audit, I*, Δδ_eq, фазы, WP-постпроцесс
python COSY/sct_study/py/run_offline_analysis.py

# 2) FOX (токи I* и контрольные точки)
python COSY/sct_study/py/generate_fox.py

# 3) COSY (нужен cosy.exe)
python COSY/src/run/run_cosy.py --pre
python COSY/sct_study/py/run_cosy_jobs.py --job validate
python COSY/sct_study/py/run_cosy_jobs.py --job working_point
python COSY/sct_study/py/run_cosy_jobs.py --job track --smoke
# полный пилот:
python COSY/sct_study/py/run_cosy_jobs.py --job track --full --stem magnetic_2

# 4) Postprocess трекинга
python COSY/sct_study/py/analyze_tracking.py
```

Для чтения уже посчитанных артефактов достаточно открыть `analysis.ipynb` — **новый трекинг не нужен**.

## Зафиксированные решения

- Только FR0 (без fringe), без FIT.
- Structures: `magnetic_2`…`magnetic_5`.
- Пилот трекинга: `magnetic_2` (лучший `cond(R_int)`).
- Multi-turn: `INJECT+TR`; анализ индивидуальных \(\overline D_i-\overline D_{\mathrm{ref}}\) и \(C(n)\). Поле `dnu_s` сохраняется только как диагностический алиасированный наклон фазы.

См. [`WORKLOG.md`](WORKLOG.md) для статусов и чисел.
