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

## Audit mapping

- **Статус:** verified (offline)
- **Команда:** `python COSY/sct_study/py/audit_mapping.py --mode both --presentation`
- **Выход:** `dat/<stem>/audit[_FR0|_FR3].json`, `plots/<stem>_grid_vs_istar*.png`

| stem | mode | cond | I* | inside box | d_min | d_min/‖ΔI‖ | plot |
|------|------|------|----|------------|-------|------------|------|
| magnetic_2 | FR0 | 34.5 | [0.004110603226596334, -0.003183980136620862, -0.004030542798297506] | False | 0.004102 | 0.452 | magnetic_2_grid_vs_istar_FR0.png |
| magnetic_2 | FR3 | 6.5 | [0.0036739260676598704, -0.002248732880017494, -0.0028758672454691765] | False | 0.005179 | 0.322 | magnetic_2_grid_vs_istar_FR3.png |
| magnetic_3 | FR0 | 35.9 | [0.027589992513895895, -0.02835137491508242, -0.02307796576882895] | False | 0.02863 | 1.54 | magnetic_3_grid_vs_istar_FR0.png |
| magnetic_3 | FR3 | 17.1 | [0.022567587446025765, -0.015348925536200068, -0.012610618398778464] | False | 0.01991 | 0.63 | magnetic_3_grid_vs_istar_FR3.png |
| electrostatic | FR0 | 47.5 | [-0.21696938660284196, 13.890952090718532, -9.517714055996425] | False | 9.907 | 2.99 | electrostatic_grid_vs_istar_FR0.png |
| electrostatic | FR3 | 47.3 | [-0.18701060835270725, 13.82789008444148, -9.574775398850074] | False | 12.36 | 3.73 | electrostatic_grid_vs_istar_FR3.png |
| Nuclotron_8 | FR0 | 250.4 | [-0.2655065733888983, 0.40924106765036267, -0.1405443623956044] | False | 0.4235 | 10.3 | Nuclotron_8_grid_vs_istar_FR0.png |
| Nuclotron_8 | FR3 | 259.1 | [-0.2888320975832907, 0.40350851394702963, -0.13231398101268288] | False | 0.4699 | 5.7 | Nuclotron_8_grid_vs_istar_FR3.png |
| Nuclotron_16 | FR0 | 406.6 | [-0.4056055937838175, 0.2640627057337158, -0.06742987492860233] | False | 0.486 | 281 | Nuclotron_16_grid_vs_istar_FR0.png |
| Nuclotron_16 | FR3 | 407.0 | [-0.39712954925097527, 0.25665487025383316, -0.06589476655127777] | False | 0.4761 | 115 | Nuclotron_16_grid_vs_istar_FR3.png |

## Audit mapping

- **Статус:** verified (offline)
- **Команда:** `python COSY/sct_study/py/audit_mapping.py --mode both --presentation`
- **Выход:** `dat/<stem>/audit[_FR0|_FR3].json`, `plots/<stem>_grid_vs_istar*.png`

| stem | mode | cond | I* | inside box | d_min | d_min/‖ΔI‖ | plot |
|------|------|------|----|------------|-------|------------|------|
| magnetic_2 | FR0 | 34.5 | [0.004110603226596334, -0.003183980136620862, -0.004030542798297506] | False | 0.004102 | 0.452 | magnetic_2_grid_vs_istar_FR0.png |
| magnetic_2 | FR3 | 6.5 | [0.0036739260676598704, -0.002248732880017494, -0.0028758672454691765] | False | 0.005179 | 0.322 | magnetic_2_grid_vs_istar_FR3.png |
| magnetic_3 | FR0 | 35.9 | [0.027589992513895895, -0.02835137491508242, -0.02307796576882895] | False | 0.02863 | 1.54 | magnetic_3_grid_vs_istar_FR0.png |
| magnetic_3 | FR3 | 17.1 | [0.022567587446025765, -0.015348925536200068, -0.012610618398778464] | False | 0.01991 | 0.63 | magnetic_3_grid_vs_istar_FR3.png |
| electrostatic | FR0 | 47.5 | [-0.21696938660284196, 13.890952090718532, -9.517714055996425] | False | 9.907 | 2.99 | electrostatic_grid_vs_istar_FR0.png |
| electrostatic | FR3 | 47.3 | [-0.18701060835270725, 13.82789008444148, -9.574775398850074] | False | 12.36 | 3.73 | electrostatic_grid_vs_istar_FR3.png |
| Nuclotron_8 | FR0 | 247.5 | [-0.2607032486637308, 0.4026997487003092, -0.13923315901443303] | False | 0.4153 | 10.1 | Nuclotron_8_grid_vs_istar_FR0.png |
| Nuclotron_8 | FR3 | 256.1 | [-0.2836773320598912, 0.39701078427973296, -0.13110463209599274] | False | 0.4255 | 5.16 | Nuclotron_8_grid_vs_istar_FR3.png |
| Nuclotron_16 | FR0 | 406.4 | [-0.4061375332269248, 0.2635254322482896, -0.06743512249165767] | False | 0.4861 | 281 | Nuclotron_16_grid_vs_istar_FR0.png |
| Nuclotron_16 | FR3 | 406.9 | [-0.39759285030783464, 0.25613362860795624, -0.06588406304131017] | False | 0.4746 | 115 | Nuclotron_16_grid_vs_istar_FR3.png |
