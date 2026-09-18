# MOCASSIN Stout Atomic Data Converter

A general, automated toolchain to ingest official [Cloudy Stout atomic database](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/StoutData) files (`.nrg`, `.tp`, `.coll`) and generate native MOCASSIN atomic collision and radiative data files (`data/*.dat`).

---

## 1. Overview & Motivation

MOCASSIN calculates thermal balance and collisionally excited line (CEL) cooling using precomputed atomic data tables located in the `data/` directory. Historically, 18 ions in MOCASSIN dated from the 1980s and 1990s (e.g. Pelan & Berrington 1995; Saraph 1986) because they were absent from CHIANTI. Furthermore, for key cooling ions (like $[\mathrm{O\ III}]$, $[\mathrm{O\ II}]$, $[\mathrm{Ne\ III}]$, $[\mathrm{S\ III}]$), the Cloudy team adopted the **Stout** database to incorporate modern $R$-matrix calculations and avoid coronal approximation limitations.

This package provides a standalone, production-grade Python converter ([`generate_stout_atomic_data.py`](file:///home/hmonteiro/software/mocassin_HM_2026/tools/stout_converter/generate_stout_atomic_data.py)) and verification suite ([`test_stout_converter.py`](file:///home/hmonteiro/software/mocassin_HM_2026/tools/stout_converter/test_stout_converter.py)) that:
1. **Bridges Cloudy Stout with MOCASSIN**: Seamlessly converts Stout's `.nrg` (energy levels), `.tp` (transition probabilities), and `.coll` (collision strengths/rates) into MOCASSIN's exact Fortran 90 format.
2. **Modernizes Legacy 1980s/1990s Ions**: Provides cutting-edge atomic data for ions missing from CHIANTI, such as:
   - **$\mathrm{F\ IV}$** (`fiv.dat`): Upgrades 1994 Lennon & Burke data to modern **ADF04 (2025-06-18)** data with complete transition probabilities.
   - **$\mathrm{Mg\ I}$** (`mgi.dat`): Upgrades 1986 Saraph JAJOM data to Barklem (2012) & Osorio (2015).
   - **$\mathrm{Cl\ IX}$** (`clix.dat`): Upgrades to Berrington, Saraph & Tully (1998).
   - **$\mathrm{Ar\ II}$, $\mathrm{Ar\ VI}$, $\mathrm{Ca\ IV}$, $\mathrm{K\ III}$, $\mathrm{K\ VII}$, $\mathrm{P\ III}$, $\mathrm{Sc\ V}$, $\mathrm{Ti\ VI}$, $\mathrm{V\ VII}$**: Modernized with NIST ASD energies and evaluated collision strengths.
3. **Implements Detailed Balance Conversion**: Automatically converts Stout `RATE ELECTRON` de-excitation rate coefficients $q_{ul}(T)$ ($\text{cm}^3\,\text{s}^{-1}$) into thermally averaged collision strengths $\Upsilon(T)$ via:
   $$\Upsilon_{lu}(T) = \frac{q_{ul}(T) \cdot g_u \sqrt{T}}{8.6291 \times 10^{-6}}$$
4. **Guarantees Fortran 90 Parser Compatibility**:
   - Enforces array bounds ($N_{\rm levs} \le 17$, matching `nForLevels = 17` in `source/constants_mod.f90`).
   - Enforces a strictly positive floor of $1.0\times 10^{-30}$ to guard against premature Fortran loop termination on $q_x == 0.0$.
   - Sums multiple radiative decay channels (e.g. $A_{\rm M1} + A_{\rm E2}$) per transition pair.
   - Resolves the Vanadium filename bug by generating both `vvii.dat` (specified in `fileNames.dat`) and compatibility alias `v7.dat`.

---

## 2. Supported Legacy Ions

| Ion | MOCASSIN File | Levels | Old Source | Stout (Cloudy c25.00) Source |
|:---|:---|:---:|:---|:---|
| **$\mathrm{Ar\ II}$** | `data/arii.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{Ar\ VI}$** | `data/arvi.dat` | 2 | Pequignot & Aldrovandi (1986) | NIST ASD / Saraph & Storey (1996) |
| **$\mathrm{Ca\ IV}$** | `data/caiv.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{Cl\ IX}$** | `data/clix.dat` | 2 | Saraph & Tully (1994) | NIST ASD / Berrington, Saraph & Tully (1998) |
| **$\mathrm{F\ II}$** | `data/fii.dat` | 5 | Butler & Zeippen (1994) | NIST ASD / Butler & Zeippen (1994) |
| **$\mathrm{F\ IV}$** | `data/fiv.dat` | 5 | Lennon & Burke (1994) (missing A) | NIST ASD / **ADF04 (2025-06-18)** |
| **$\mathrm{K\ III}$** | `data/kiii.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{K\ VII}$** | `data/kvii.dat` | 2 | Saraph & Storey (1996) | NIST ASD / Saraph & Storey (1996) |
| **$\mathrm{Mg\ I}$** | `data/mgi.dat` | 5 | Saraph (1986) / Mendoza (1983) | NIST ASD / **Barklem (2012) & Osorio (2015)** |
| **$\mathrm{P\ III}$** | `data/piii.dat` | 2 | Saraph & Storey | NIST ASD / Krueger & Czyzak (1970) |
| **$\mathrm{Sc\ V}$** | `data/scv.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{Ti\ VI}$** | `data/tivi.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{V\ VII}$** | `data/vvii.dat` / `v7.dat` | 2 | Pelan & Berrington (1995) | NIST ASD / Pelan & Berrington (1995) |
| **$\mathrm{S\ III}$** | `data/siii.dat` | 5 | Mendoza (1983) | NIST ASD / Hudson et al. (2012) |

---

## 3. CLI Usage & Common Workflows

```bash
# General syntax
/home/hmonteiro/miniforge3/bin/python tools/stout_converter/generate_stout_atomic_data.py [OPTIONS]
```

### Options

| Flag | Type | Default | Description |
|:---|:---:|:---:|:---|
| `--stout-dir` | Path | Auto | Path to Cloudy Stout directory (auto-detects local Cloudy c25.00) |
| `--ions`, `--species` | String | `legacy` | Comma-separated ions (e.g. `f_4,mg_1`), `legacy`, or `all` |
| `--output-dir` | Path | `data.stout` | Directory where `.dat` files are written |
| `--max-levels` | Int | `17` | Array bounds cap ($N_{\rm levs} \le 17$) |
| `--nlevs` | String | None | Per-ion level overrides (e.g. `f_4=5,o_3=15`) |
| `--resample-grid` | Flag | `False` | Resample collision strengths onto standard 23-point temperature grid |
| `--copy-to-data` | Flag | `False` | Copy generated files directly into `data/` |

### Examples

#### Example 1: Convert all 14 Legacy Ions
```bash
/home/hmonteiro/miniforge3/bin/python tools/stout_converter/generate_stout_atomic_data.py \
    --ions legacy \
    --output-dir data.stout
```

#### Example 2: Convert Primary Nebular Coolers
```bash
/home/hmonteiro/miniforge3/bin/python tools/stout_converter/generate_stout_atomic_data.py \
    --ions o_3,o_2,ne_3,s_2,c_2 \
    --output-dir data.stout
```

#### Example 3: Deploy Directly into Production `data/`
```bash
/home/hmonteiro/miniforge3/bin/python tools/stout_converter/generate_stout_atomic_data.py \
    --ions legacy \
    --copy-to-data
```

---

## 4. Verification Suite

Run the automated regression test suite:

```bash
/home/hmonteiro/miniforge3/bin/python tools/stout_converter/test_stout_converter.py
```

The test suite validates:
1. Mathematical exactness of `RATE ELECTRON` to $\Upsilon(T)$ detailed balance inversion.
2. Bit-for-bit consistency with the verified production $[\mathrm{S\ III}]$ benchmark.
3. Standalone `gfortran` runtime execution: all 14 legacy files parse through the MOCASSIN Fortran reader with 0 errors.
4. Prevention of premature loop termination ($q_x \ge 1.0\times 10^{-30}$).
5. Compliance with array bounds ($N_{\rm levs} \le 17$).
6. Grid resampling to standard 23-point temperature grid.
