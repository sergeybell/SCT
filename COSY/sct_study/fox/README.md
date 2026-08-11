# `fox/` — входные файлы COSY

## Шаблоны

| Файл | Роль |
|------|------|
| `_template_validate.fox` | chrom, η₁, \(\Delta\nu_s\) при заданных токах |
| `_template_working_point.fox` | \(Q_x,Q_y,\nu_s,\gamma G\) → `wp_<tag>.dat` |
| `_template_track.fox` | `INJECT` + `TR` → PRAY / TRPRAY / TRPSPI |

Плейсхолдеры заполняет `py/generate_fox.py`.

## Сгенерированные файлы

- `validate_<stem>_{natural,Istar}.fox`
- `working_point_<stem>_<tag>.fox`
- `track_<stem>_<tag>_{smoke,full,dense}.fox` (пилот: `magnetic_2`)

## Ансамбль `INJECT`

1. Ray 0 — служебный COSY; ray 1 — reference (\(X=\ldots=D=0\)).
2. Группы X / Y / D по `num_per_group` (linspace; в середине группы амплитуда 0 — классификация по порядку inject, не по |X|).
3. Спин: `psi_deg=0` → \(S=(0,0,1)\) (горизонтальная/продольная фаза `atan2(Sx,Sz)`).
4. Сохранение: `TR NTURN SAVE_EVERY ...` с явным `SAVE_EVERY` (dense: 2).
