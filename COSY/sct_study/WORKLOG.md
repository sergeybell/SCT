# WORKLOG — SCT study

## Scope
- Structures: `magnetic_2`…`magnetic_5`
- Mode: FR0 only (no fringe), no FIT
- Folder: `COSY/sct_study/`
- Pilot tracking: `magnetic_2`

---

## 1. Setup — verified
- Created `README.md`, `WORKLOG.md`, `manifest.json`, `config/`, `fox/`, `py/`, `dat/`, `plots/`

## 2. Audit mapping FR0 — verified
Command: `python COSY/sct_study/py/audit_mapping.py`

| stem | cond(R_int) | I* inside box? | min\|I*-grid\| |
|------|-------------|----------------|----------------|
| magnetic_2 | 34.5 | no (SGx2<0) | 0.0041 |
| magnetic_3 | 35.9 | no | 0.029 |
| magnetic_4 | 180.6 | no, far | 0.34 |
| magnetic_5 | 69.7 | no, far | 0.24 |

- Linearity R²≈1 for chrom and spin-link.
- Outputs: `dat/<stem>/audit.json`, `plots/<stem>_grid_vs_istar.png`

## 3. Coefficient conventions — verified
Command: `python COSY/sct_study/py/validate_coefficients.py`
- γ=1.143914, G=−0.142987, γG≈−0.163565
- η₁ rebuild from α₀,α₁ matches integrals.dat
- WARNING: MCM2 factorial differs (`mapping.fox` vs `coherence_opt.fox`)
- WARNING: `fACCLEN=141` ≠ geometric L (117–134 m)
- Output: `dat/coefficient_conventions.json`

## 4. Zero point I* — verified (COSY)
Commands:
- `python COSY/sct_study/py/analyze_zero_point.py`
- `python COSY/sct_study/py/run_cosy_jobs.py --job validate`

magnetic_2 direct COSY @I*:
- SGx1=0.0041106, SGx2=−0.0031840, SGy1=−0.0040305
- chrom_x≈−9.3e−13, chrom_y≈−8.1e−13, eta1≈−1.3e−13
- Model chrom error ~1e−12 → **ξ and η₁ действительно обнуляются**

Same validation run for magnetic_3/4/5.

## 5. Working points — verified
Commands:
- `python COSY/sct_study/py/run_cosy_jobs.py --job working_point`
- `python COSY/sct_study/py/analyze_working_point.py`

magnetic_2 @I*: Qx≈0.44628, Qy≈0.47668, ν_s≈0.163565, γG≈−0.163565  
Nearest resonances far (min distance ~0.28).  
Plot: `plots/working_points_summary.png`

## 6. Phase diagrams — verified (μ proxy)
Command: `python COSY/sct_study/py/analyze_phase_advance.py`
- 2D polar: length=|K₂|, angle=μ=2πQ·s/C (proxy; real Q from WP)
- Plots: `plots/<stem>_sext_phase_{x,y}.png`

## 7. Δδ_eq panels — verified
Command: `python COSY/sct_study/py/analyze_delta_eq.py`
- Theory Δδ_eq crosses 0 at ξ=0 by construction
- Plots: `plots/<stem>_delta_eq_panels.png`, `plots/<stem>_map_dnu_panels.png`
- Summary: `plots/summary_delta_eq_suppression.png`

## 8. Multi-turn INJECT+TR — verified (pilot magnetic_2)
Commands:
- `python COSY/sct_study/py/generate_fox.py`
- `python COSY/sct_study/py/run_cosy_jobs.py --job track --smoke|--full --stem magnetic_2`
- `python COSY/sct_study/py/analyze_tracking.py`

~29 rays (ref + X/Y/D groups), smoke=200 / full=2000 turns.

Per-group RMS Δδ_eq (full):

| tag | RMS_X | RMS_Y | RMS_D | C_final |
|-----|-------|-------|-------|---------|
| natural | 6.5e−9 | 3.1e−8 | 1.5e−5 | 9.2e−4 |
| **Istar** | **5.7e−13** | **7.6e−12** | 1.5e−5 | 4.6e−4 |
| ctrl_xi_x | 6.9e−8 | 3.2e−7 | 1.5e−5 | 1.2e−2 |

**Physics conclusion:**
1. Mapping + I* correctly zeros ξ_x, ξ_y, η₁ (COSY confirmed).
2. Multi-turn path-lengthening proxy for X/Y groups drops by ~4 orders at I*.
3. D-group unchanged (synchrotron motion).
4. Spin coherence C(n) on this diagnostic ensemble does **not** improve at I* — next checks: fix fACCLEN to geometric L, spin-plane definition, matched beam ensemble.

Plots: `plots/magnetic_2_track_*.png`, `plots/magnetic_2_track_compare_full.png`

---

## Reproduce
```bash
python COSY/sct_study/py/run_offline_analysis.py
python COSY/src/run/run_cosy.py --pre
python COSY/sct_study/py/run_cosy_jobs.py --job validate
python COSY/sct_study/py/run_cosy_jobs.py --job working_point
python COSY/sct_study/py/run_cosy_jobs.py --job track --full --stem magnetic_2
python COSY/sct_study/py/analyze_working_point.py
python COSY/sct_study/py/analyze_tracking.py
```

---

## 9. Documentation canvas — verified
- README в `config/`, `py/`, `fox/`, `dat/` (+ stems), `plots/`; расширен корневой `README.md` (словарь smoke/full, C(n)).
- Единый notebook: `analysis.ipynb` (физика + воспроизведение без нового COSY track).
- Полярные диаграммы: свой масштаб на stem + результирующий вектор Σ; поле `resultant` в `sext_phase.json`.

## Phase diagrams

- **Статус:** prepared (s/C·Q proxy for μ; update after working_point COSY run)
- **Команда:** `python COSY/sct_study/py/analyze_phase_advance.py`
- Диаграмма: длина = |K₂|, угол = μ; знак K₂ → +π; чёрный пунктир = результирующий Σ.
- Масштаб радиуса — свой для каждой структуры.
- magnetic_2: 10 magnets; |Σ|_x=0.01527, |Σ|_y=0.01457; magnetic_2_sext_phase_x.png, magnetic_2_sext_phase_y.png
- magnetic_3: 15 magnets; |Σ|_x=0.08451, |Σ|_y=0.09117; magnetic_3_sext_phase_x.png, magnetic_3_sext_phase_y.png
- magnetic_4: 20 magnets; |Σ|_x=0.2824, |Σ|_y=0.3531; magnetic_4_sext_phase_x.png, magnetic_4_sext_phase_y.png
- magnetic_5: 25 magnets; |Σ|_x=1.302, |Σ|_y=1.316; magnetic_5_sext_phase_x.png, magnetic_5_sext_phase_y.png
