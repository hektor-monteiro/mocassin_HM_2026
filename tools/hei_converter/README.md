# MOCASSIN He I Recombination Line Converter

This tool generates modern Case B He I recombination line coefficients for MOCASSIN based on the comprehensive calculations by **Porter, Ferland, Storey, & Detisch (2012, MNRAS, 425, L28; erratum 2013, MNRAS, 433, L89)**.

---

## 1. Physics Background

In MOCASSIN, the volume emissivity of bound-bound He I recombination lines is represented by the parametric form:
$$4\pi j(\lambda) = \left[ A \cdot T_4^b \cdot \exp(c / T_4) \right] \cdot N_e N(\mathrm{He}^+)$$
where $T_4 = T_e / 10^4\ \mathrm{K}$.

- **Legacy Source**: Benjamin, Skillman, & Smits (1999, ApJ, 514, 307), based on Smits (1996).
- **Modern Source**: Porter et al. (2012, 2013), incorporating $n \le 100$ resolved energy levels, updated radiative transitions, and ab initio $R$-matrix electron-impact collision strengths (Bray et al. 2000; Ballance et al. 2006).

### Improvements in Porter et al. (2012/2013)
1. **Collisional Transfers**: Accurate inclusion of $l$-mixing and spin-changing collisions between $n=3$ terms ($3^3D \leftrightarrow 3^1D$).
2. **High-Density Enhancements**: Lines such as $[\mathrm{He\ I}]\ 6678\ \text{Å}$ and $7065\ \text{Å}$ receive up to $+30\%\text{--}70\%$ enhancements at densities $N_e \ge 10^4\ \mathrm{cm}^{-3}$.
3. **Low-Density Consistency**: In low-density photoionized environments ($N_e \sim 10^2\ \mathrm{cm}^{-3}$), optical lines match the legacy dataset to within $\le 1\%$, guaranteeing smooth continuity with historical benchmarks.

---

## 2. Directory Structure

- `generate_porter_hei_data.py`: Core fitting engine and MOCASSIN data file generator.
- `test_hei_converter.py`: Verification test suite (extraction, non-linear fitting, Fortran reader compatibility, physical accuracy).
- `README.md`: Documentation and usage instructions.

---

## 3. Usage

To generate `data/HeIrecLines.dat`:

```bash
/home/hmonteiro/miniforge3/bin/python tools/hei_converter/generate_porter_hei_data.py \
    --output data/HeIrecLines.dat \
    --verbose
```

To run the verification test suite:

```bash
/home/hmonteiro/miniforge3/bin/python tools/hei_converter/test_hei_converter.py
```
