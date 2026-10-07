# MOCASSIN - PyNeb Emission Line Intensity & ADF Calculator (`mocassin_pyneb.py`)

A Python tool to compute 3D integrated emission line intensities from MOCASSIN photoionization models using **PyNeb** and the **CHIANTI** atomic database.

---

## Key Capabilities & Architecture

### 1. CHIANTI v10 Integration with Exact MOCASSIN Level Counts (`NLEVS`)
- **Database**: Configured to use the official CHIANTI v10.1 database (extracted at `/home/hmonteiro/software/chianti_v10.1` or passed via `--chianti-dir`).
- **Energy Level Matching**: Dynamically reads MOCASSIN's exact level counts (`NLEVS`) from `data/fileNames.dat` (or uses built-in default mapping) across all 116 ions.
  - Examples: [O III] (15 levels), [O II] (13 levels), [O I] (7 levels), [N II] (15 levels), [N III] (10 levels), [C II] (10 levels), [C III] (10 levels), [S II] (5 levels), [S III] (5 levels), [Ne III] (9 levels), [Fe II] (142 levels).
  - This ensures PyNeb calculates collisional line emissivities using the same multi-level statistical equilibrium matrix as MOCASSIN.

### 2. Collisional (CEL) and Recombination (ORL) Line Support
- **Collisional Excitation Lines (CELs)**:
  - Evaluated using `pyneb.Atom(elem, spec, NLevels=NLEVS)` with CHIANTI v10 atomic data.
  - Volume integrated over local electron density $N_e$, electron temperature $T_e$, emitting ion density $N(X^{i+})$, and cell volume element $dV$:
    $$L_{\text{CEL}} = \sum_{\text{cells}} \epsilon(T_e, N_e) \cdot N_e \cdot N(X^{i+}) \cdot dV \quad [\text{erg s}^{-1}]$$
- **Optical Recombination Lines (ORLs)**:
  - Evaluated using `pyneb.RecAtom(elem, spec)`.
  - Volume integrated over the **parent recombining ion density** $N(X^{(i+1)+})$:
    $$L_{\text{ORL}} = \sum_{\text{cells}} \epsilon_{\text{eff}}(T_e, N_e) \cdot N_e \cdot N(X^{(i+1)+}) \cdot dV \quad [\text{erg s}^{-1}]$$
  - Examples: O II ORLs integrate over $N(O^{2+})$, C II ORLs over $N(C^{2+})$, N II ORLs over $N(N^{2+})$.

### 3. 3D Grid & Abundance Reconstruction
- Reconstructs Cartesian coordinates $(x, y, z)$ and differential volume elements $dV = dx \cdot dy \cdot dz$ [$\text{cm}^3$] from MOCASSIN outputs:
  - `grid0.out` (cell coordinates and spatial grid structure)
  - `grid1.out` (electron temperatures $T_e$ and electron densities $N_e$)
  - `grid2.out` (total gas density and elemental abundances, supporting `multiChemistry` spatial variations)
  - `grid3.out` (ionic fraction tables for each grid cell)
- Automatically takes into account octant reflection symmetry (`lgSymmetricXYZ` factor $2.0\times$ in half-plane or octant models).

### 4. Planetary Nebula ADF Diagnostics
- Automatically computes ionic abundances relative to $H^+$:
  - $(X^{2+}/H^+)_{\text{CEL}}$ from dominant optical/far-IR collisionally excited lines.
  - $(X^{2+}/H^+)_{\text{ORL}}$ from optical recombination lines.
- Computes the **Abundance Discrepancy Factor**:
  $$\text{ADF}(X^{2+}) = \frac{(X^{2+}/H^+)_{\text{ORL}}}{(X^{2+}/H^+)_{\text{CEL}}}$$
  for $O^{2+}$ ([O III] 5007 vs O II 4649 / multiplet 1), $C^{2+}$ (C III] 1908 vs C II 4267), and $N^{2+}$ ([N III] 57 $\mu$m vs N II 5680).

### 5. Direct Comparison Against MOCASSIN `lineFlux.out`
- Cross-matches lines directly against MOCASSIN's analytical and formal solutions in `output/lineFlux.out`.
- Side-by-side terminal report and CSV export reporting percentage residuals:
  $$\Delta(\%) = 100 \times \frac{I_{\text{PyNeb}} - I_{\text{MOCASSIN}}}{I_{\text{MOCASSIN}}}$$

---

## Test Emission Line List

The test line list is provided at `benchmarks/test_Problem/input/lines_plot_and_adf.dat`:

1. **All 47 Lines from `plot.in`**:
   - **H & He**: H$\alpha$ (6562.8 Å), H$\beta$ (4861.3 Å), H$\gamma$ (4340.5 Å), He I (5876, 6678, 4471 Å), He II (1640, 4686, 5411 Å).
   - **Carbon**: [C I] 8727 Å, [C II] 158 $\mu$m, C III] 1908/1906 Å.
   - **Nitrogen**: [N I] 5198/5200 Å, [N II] 6548/6584/5755 Å, [N II] 122 $\mu$m, [N III] 57 $\mu$m.
   - **Oxygen**: [O I] 6300/5577 Å, [O I] 63/145 $\mu$m, [O II] 3726/3729/7320/7330 Å, [O III] 5007/4363 Å, [O III] 52/88 $\mu$m.
   - **Neon & Sulfur**: [Ne III] 3869/3967 Å, [S II] 6716/6731 Å, [S III] 6312 Å, [S III] 18.7 $\mu$m, [S III] 9069 Å.
   - **Chlorine & Argon**: [Cl III] 5518/5538 Å, [Ar III] 7136 Å, [Ar IV] 4711/4740 Å.

2. **Key Planetary Nebula ADF Optical Recombination Lines (ORLs)**:
   - **O II Multiplet 1 (V1)** ($3s\,^4\text{P} - 3p\,^4\text{D}^\circ$): 4649.13 Å (principal), 4650.84 Å, 4661.63 Å, 4676.23 Å, 4673.73 Å, 4641.81 Å, 4638.86 Å.
   - **Other O II ORLs**: 4089.29 Å ($3d\,^4\text{F} - 4f\,\text{G}[5]^\circ$, M48), 4596.18 Å ($3s'\,^2\text{D} - 3p'\,^2\text{F}^\circ$, M15).
   - **C II ORLs**: 4267.00 Å ($3d\,^2\text{D} - 4f\,^2\text{F}^\circ$, classic ADF line), 6578.05 Å ($3s\,^2\text{S} - 3p\,^2\text{P}^\circ$).
   - **N II ORLs**: Multiplet 3 (5679.56 Å, 5666.63 Å, 5676.02 Å), 4630.54 Å (M5), 4041.31 Å (M39).

3. **Permitted Recombination Lines Contaminating Auroral / Diagnostic Lines**:
   - **O II Multiplet 2 (V2)** ($3s\,^4\text{P} - 3p\,^4\text{P}^\circ$):
     - **O II 4366.89 Å** ($5/2 - 5/2$) & **4366.53 Å** ($3d\,^4\text{F} - 4f\,\text{D}[3]^\circ$): Directly blend with and contaminate the weak auroral [O III] $\lambda 4363.21$ line in low/intermediate resolution spectra ($R \lesssim 1200$), artificially boosting the measured flux and causing overestimates in derived $T_e([\text{O III}])$.
     - **O II 4349.43 Å** ($3/2 - 5/2$): The clean reference line of Multiplet 2, used observationally to deblend and subtract the 4366.89 Å contribution.
     - **Other V2 lines**: 4317.14 Å, 4319.63 Å, 4336.86 Å, 4345.56 Å, 4325.76 Å, 4357.25 Å.
   - **Other nearby permitted lines**: O II 4369.27 Å (M65), O II 4371.61 Å (M44), N II 5686.21/5710.77 Å (near [N II] 5755), and O II 7318.99 Å (blends with [O II] 7320).

---

## Verification on `benchmarks/test_Problem/`

### Execution Command
```bash
/home/hmonteiro/miniforge3/envs/astro/bin/python accessories/mocassin_pyneb.py \
  -m benchmarks/test_Problem \
  -l benchmarks/test_Problem/input/lines_plot_and_adf.dat \
  -o benchmarks/test_Problem/output/pyneb_intensities.csv
```

### 1. Comparison with MOCASSIN Solutions
- **Total $H\beta$ Luminosity**:
  - PyNeb: $3.30296 \times 10^{35} \text{ erg s}^{-1}$ (Storey & Hummer 1995 Case B)
  - MOCASSIN: $3.20255 \times 10^{35} \text{ erg s}^{-1}$ (Brocklehurst / Hummer power-law approximation)
  - Denominator offset: $+3.14\%$
- **Absolute Line Luminosities**:
  - Absolute luminosities (e.g. $[O\text{ III}]\,5007$: PyNeb $= 5.6193 \times 10^{36} \text{ erg s}^{-1}$ vs MOCASSIN $= 5.6169 \times 10^{36} \text{ erg s}^{-1}$) agree to **$< 0.05\%$**.
- **Line Ratios ($I / I(H\beta)$)**:
  - H$\alpha$ / H$\beta$: PyNeb = 2.798 vs MOCASSIN = 2.798 ($\Delta = +0.02\%$)
  - H$\gamma$ / H$\beta$: PyNeb = 0.472 vs MOCASSIN = 0.472 ($\Delta = -0.02\%$)
  - [O III] 5007 / H$\beta$: PyNeb = 17.01 vs MOCASSIN = 17.54 ($\Delta = -3.00\%$)
  - [O III] 4363 / H$\beta$: PyNeb = 0.288 vs MOCASSIN = 0.297 ($\Delta = -3.03\%$)
  - [N II] 6584 / H$\beta$: PyNeb = 2.311 vs MOCASSIN = 2.383 ($\Delta = -3.04\%$)
  - [Ne III] 3869 / H$\beta$: PyNeb = 2.406 vs MOCASSIN = 2.481 ($\Delta = -3.00\%$)
  - [S II] 6731 / H$\beta$: PyNeb = 0.661 vs MOCASSIN = 0.682 ($\Delta = -3.00\%$)

### 2. ADF Planetary Nebula Diagnostics
```text
===========================================================================
Abundance Discrepancy Factor (ADF) Planetary Nebula Diagnostics
===========================================================================
Parent Ion      CEL Line        ORL Line         (X2+/H+) CEL   (X2+/H+) ORL   ADF     
---------------------------------------------------------------------------
O(2+) / H(+)    [O III] 5007    O II 4649.13          2.161e-04      2.178e-04     1.01
  -> Total O II Multiplet 1 (V1) Sum I(V1)/I(Hb) = 2.7984e-03 (0.2798 on Hb=100)
C(2+) / H(+)    C III] 1908     C II 4267.0           1.432e-04      1.449e-04     1.01
N(2+) / H(+)    [N III] 57m     N II 5679.56          5.405e-05      3.385e-05     0.63
===========================================================================
```
*Because `test_Problem` is a chemically homogeneous photoionization model, the theoretical ADF for oxygen and carbon is expected to be ~1.0. PyNeb yields $\text{ADF} = 1.01$, validating atomic consistency across the 3D volume.*

### 3. Auroral Line Recombination Contamination Assessment
```text
==========================================================================================
Auroral Line Contamination Assessment: [O III] 4363 by Permitted O II ORLs
==========================================================================================
  Forbidden [O III] 4363.21 Intensity (Hb=100)           :    28.7904
  Permitted O II 4366.89 (Multiplet 2, 5/2-5/2) (Hb=100)  :     0.0151
  Permitted O II 4366.53 (3d 4F - 4f D[3]0) (Hb=100)     :     0.0011
  Total Permitted O II Blend near 4363 A (Hb=100)        :     0.0162
  Total O II Multiplet 2 (V2) Sum (Hb=100)               :     0.0916
  Blend Contamination Fraction                           :     0.056% of [O III] 4363
  Uncontaminated Ratio R([O III]) = (4959+5007)/4363     :      79.53 -> Te = 13884.3 K
  Contaminated Ratio (including unresolved O II blend)   :      79.49 -> Te = 13884.3 K (delta Te =  +0.0 K)
==========================================================================================
```
*In low-temperature ($T_e \lesssim 9,000\text{ K}$) nebulae or objects with large Abundance Discrepancy Factors where ORLs are heavily enhanced and [O III] 4363 is weak, this unresolved blend can contribute several percent or tens of percent to the measured 4363 flux, leading to significant overestimates in the derived [O III] electron temperature.*

---

## Command-Line Options Reference

```text
usage: mocassin_pyneb.py [-h] --model-dir MODEL_DIR [--line-list LINE_LIST]
                         [--plot-in PLOT_IN] [--chianti-dir CHIANTI_DIR]
                         [--mocassin-dir MOCASSIN_DIR] [--output OUTPUT]
                         [--adf] [--compare-mocassin] [--no-compare-mocassin]
```

| Option | Flag | Description |
|---|---|---|
| `--model-dir` | `-m` | Path to MOCASSIN model directory containing `output/` and `input/` (**Required**). |
| `--line-list` | `-l` | Path to emission line list (`Ion Wavelength [Type] [Label]`). Defaults to `lines_plot_and_adf.dat` or `plot.in`. |
| `--plot-in` | `-p` | Path to MOCASSIN `plot.in` format line list. |
| `--chianti-dir`| `-c` | Path to official CHIANTI database directory. Default searches `/home/hmonteiro/software/chianti_v10.1` or standard env. |
| `--mocassin-dir` | | Path to MOCASSIN source directory to read level counts `NLEVS` from `data/fileNames.dat`. |
| `--output` | `-o` | Output file path to save CSV results. |
| `--adf` | | Display Abundance Discrepancy Factor diagnostics for PNe (Default: `True`). |
| `--compare-mocassin` | | Compare against MOCASSIN `output/lineFlux.out` (Default: enabled if file exists). |
| `--no-compare-mocassin` | | Disable comparison against `lineFlux.out`. |
