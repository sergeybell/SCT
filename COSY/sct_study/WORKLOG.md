# WORKLOG — SCT study

## Scope
- Structures: `magnetic_2`…`magnetic_5`
- Mode: FR0 only (no fringe), no FIT
- Folder: `COSY/sct_study/`
- Pilot tracking: `magnetic_2`
- Spin sampling: dense \(\Delta n=2\), horizontal spin (`psi_deg=0`)

---

## 1. Mapping FR0 — verified
| stem | cond(R_int) | inside FR0 | d_min | d_min/‖ΔI‖ |
|------|-------------|------------|-------|------------|
| magnetic_2 | 34.5 | no (SGx2<0) | 0.0041 | **0.45** |
| magnetic_3 | 35.9 | no | 0.029 | 1.54 |
| magnetic_4 | 180.6 | no, far | 0.34 | 10.7 |
| magnetic_5 | 69.7 | no, far | 0.24 | 5.16 |

Linearity R²≈1. Сравнивать stems по **нормированному** расстоянию.

## 2. Coefficients — verified
- γ=1.143914, G=−0.142987, γG≈−0.163565
- η₁ rebuild matches integrals.dat
- MCM2: study uses `mapping.fox` convention
- WARNING: `fACCLEN=141` ≠ geometric L (117–134 m)

## 3. I* validate (COSY) — verified
magnetic_2 @I*: chrom_x, chrom_y, η₁ ~ 1e−12 → обнуление подтверждено.

## 4. Working points — verified
magnetic_2 @I*: Qx≈0.44628, Qy≈0.47668, ν_s≈0.163565

## 5. Twiss phase at sextupole instances — verified
- ψ=∫ds/β from existing BETAX/BETAY (linear β independent of sextupoles on closed orbit)
- rms |ψ − 2πQs/C| ~ O(1) rad → геометрический proxy **недостаточен**
- Outputs: `dat/<stem>/sext_phase_twiss.json`, `plots/<stem>_sext_phase_twiss_{x,y}.png`

## 6. Resonance scan — verified (geometric + tracking proxies)
Kinds: imperfection / intrinsic / synchrotron / combined.
- Deuteron |Gγ|≪ proton at same γ — but coherence still needs explicit check.
- Qs currently from config estimate (0.01) until RF one-turn extraction.
- Plot: `plots/resonance_scan_summary.png`

## 7. Multi-turn tracking — verified

### Nomenclature
- `mean_D_offset` = ⟨D⟩−⟨D⟩ref (not theoretical Δδ_eq)
- Spin phase requires horizontal initial spin (`psi_deg=0`)
- Sparse `full` (Δn=100) / `smoke` (Δn=10): **aliased** for spin; orbit proxy OK

### Dense (Δn=2, NTURN=2000, psi=0), magnetic_2

| tag | RMS_X | RMS_Y | RMS_D | C_final | dnu_rms | ref ν_s |
|-----|-------|-------|-------|---------|---------|---------|
| natural | 6.1e−9 | 2.9e−8 | 1.3e−5 | 1.0 | ~4e−8 | 0.16356483 |
| **Istar** | **4.4e−13** | **6.7e−12** | 1.3e−5 | 1.0 | ~4e−8 | 0.16356483 |
| ctrl_xi_x | 6.4e−8 | 3.0e−7 | 1.3e−5 | 1.0 | ~4e−8 | 0.16356483 |

**Physics:**
1. I* zeros ξ, η₁ and suppresses X/Y `mean_D_offset` by ~4 orders.
2. Dense ref spin tune matches |γG| to machine precision.
3. On this diagnostic ensemble (2000 turns, amp 1 mm / δ~1e−4) relative Δν and C(n) remain ~ideal for all tags — decoherence not resolved yet; next: longer run / matched larger emittance / fix fACCLEN.

## Reproduce
```bash
python COSY/sct_study/py/run_offline_analysis.py
python COSY/sct_study/py/generate_fox.py
python COSY/sct_study/py/run_cosy_jobs.py --job track --dense --stem magnetic_2
python COSY/sct_study/py/analyze_tracking.py
python COSY/sct_study/py/analyze_resonances.py
```

## Resonance scan

- **Статус:** verified (geometric distances + tracking proxies when present)
- **Команда:** `python COSY/sct_study/py/analyze_resonances.py`
- Plot: `resonance_scan_summary.png`
- Qs from config estimate until RF one-turn extraction is added.
- magnetic_2@Istar: imp=0.1636 (ν_s=5), intr=0.2827 (ν_s=4+Q_x), sync=0.1436, comb=0.2627
- magnetic_3@Istar: imp=0.1636 (ν_s=0), intr=0.1269 (ν_s=3-Q_x), sync=0.1436, comb=0.1069
- magnetic_4@Istar: imp=0.1636 (ν_s=5), intr=0.07527 (ν_s=6-Q_y), sync=0.1436, comb=0.05527
- magnetic_5@Istar: imp=0.1636 (ν_s=3), intr=0.06388 (ν_s=1+Q_y), sync=0.1436, comb=0.04388

## Δδ_eq + map-Δν_s panels

- **Статус:** verified (offline)
- **Команда:** `python COSY/sct_study/py/analyze_delta_eq.py`
- Сценарии ε/δ заданы в `config/study_config.json`.
- При линейной модели Δδ_eq(ξ) проходит через 0 при ξ=0 **по построению**.
- magnetic_2: magnetic_2_delta_eq_panels.png, magnetic_2_map_dnu_panels.png; Δδ_eq@I*={'X_only': -1.1131338950049443e-23, 'X_only_at_natural': -9.228135434107606e-08, 'Y_only': -1.1131338950049443e-23, 'Y_only_at_natural': -9.262616120172013e-08, 'D_only': -2.079659993557127e-24, 'D_only_at_natural': -9.532609864673429e-09}
- magnetic_3: magnetic_3_delta_eq_panels.png, magnetic_3_map_dnu_panels.png; Δδ_eq@I*={'X_only': -1.7204862495743977e-23, 'X_only_at_natural': -1.447261149502933e-07, 'Y_only': -1.7204862495743977e-23, 'Y_only_at_natural': -1.4510670504559358e-07, 'D_only': -0.0, 'D_only_at_natural': -7.739409948872235e-09}
- magnetic_4: magnetic_4_delta_eq_panels.png, magnetic_4_map_dnu_panels.png; Δδ_eq@I*={'X_only': 6.170243037032215e-23, 'X_only_at_natural': -1.9065534057985742e-07, 'Y_only': -1.234048607406443e-22, 'Y_only_at_natural': -1.9051406274606084e-07, 'D_only': 3.147101757126394e-24, 'D_only_at_natural': -6.882172082509207e-09}
- magnetic_5: magnetic_5_delta_eq_panels.png, magnetic_5_map_dnu_panels.png; Δδ_eq@I*={'X_only': -2.8745982597068286e-23, 'X_only_at_natural': -2.3470501642833182e-07, 'Y_only': -5.749196519413657e-23, 'Y_only_at_natural': -2.3763211625023798e-07, 'D_only': -7.639689592633897e-25, 'D_only_at_natural': -6.44177135024405e-09}

## OptiM Nu vs Twiss vs Q_full·s/C

- **Статус:** offline from `sext_phase_magnetic_2.txt` + I* + existing Twiss JSON
- **Команда:** `python COSY/sct_study/py/analyze_optim_sext_phase.py`
- A: |I*| with SF1>0, SF2>0, SD<0. φ=2π Nu (Nu in turns).
- magnetic_2: 10 instances; rms|ΔNu|_geom=0.008911/0.009622; rms|ΔNu|_twiss=0.2904297853494206/0.22704263963720983
