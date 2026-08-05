# `fox/` — входные файлы COSY

## Шаблоны (редактировать вручную)

| Файл | Роль |
|------|------|
| `_template_validate.fox` | Одноточечная проверка chrom, η₁, \(\Delta\nu_s\) при заданных токах |
| `_template_working_point.fox` | Рабочая точка: \(Q_x,Q_y,\nu_s,\gamma G\) → `wp_<tag>.dat` |
| `_template_track.fox` | `INJECT` + multi-turn `TR` → PRAY / TRPRAY / TRPSPI |

Плейсхолдеры вида `{{STEM}}`, `{{SGX1}}`, `{{NTURN}}` заполняет `py/generate_fox.py`.

## Сгенерированные файлы (не править руками)

Имена:

- `validate_<stem>_<tag>.fox` — `tag ∈ {natural, Istar}`
- `working_point_<stem>_<tag>.fox` — `tag ∈ {natural, Istar, ctrl_xi_x, ctrl_xi_y, ctrl_eta1}`
- `track_<stem>_<tag>_<kind>.fox` — `kind ∈ {smoke, full}` (пилот: `magnetic_2`)

## Ансамбль `INJECT` (трекинг)

В `_template_track.fox`:

1. Reference: \(X=A=Y=B=T=D=0\), спин \(S=(0,\sin\psi,\cos\psi)\), \(\psi=90^\circ\) → \(S=(0,1,0)\).
2. Группа X: `NUM` лучей с \(X\in[-x_{amp},x_{amp}]\).
3. Группа Y: то же по \(Y\).
4. Группа D: \(D = X\cdot d_{scale}\) (масштаб импульсного разброса).

Далее `TR NTURN NINT(NTURN*save_every_frac) ...` пишет орбиту и спин.
