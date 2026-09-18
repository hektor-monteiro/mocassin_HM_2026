# MOCASSIN CHIANTI Atomic Data Converter

A robust, automated toolchain to ingest official [CHIANTI atomic database](https://www.chiantidatabase.org/) releases (v7 through v11+) and generate native MOCASSIN atomic collision and radiative data files (`data/*.dat`).

---

## 1. Overview & Motivation

MOCASSIN (*MOnte CArlo SimulationS of Ionised Nebulae*) calculates thermal balance and collisionally excited line (CEL) cooling using precomputed atomic data tables located in the `data/` directory. Historically, these files were manually assembled or generated using legacy IDL scripts that became out of sync with modern CHIANTI data formats.

This package provides a standalone, production-grade Python converter ([`generate_chianti_atomic_data.py`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/tools/chianti_converter/generate_chianti_atomic_data.py)) and test suite ([`test_converter.py`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/tools/chianti_converter/test_converter.py)) that:
1. **Eliminates manual file editing**: Automatically reads CHIANTI database releases (unpacked folders, `.tar.gz` distribution archives, or direct online downloads).
2. **Reverse-engineers MOCASSIN's exact Fortran parser**: Strictly adheres to the Fortran 90 I/O specifications in [`source/hydro_mod.f90`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/source/hydro_mod.f90#L521-L663).
3. **Resolves critical Fortran parsing pitfalls**: Prevents premature loop termination caused by zero-valued collision strengths, enforces array bounds (`nForLevels = 17`), and preserves the specialized 142-level Fe II model.
4. **Expands ion coverage**: In addition to updating the 96 historical MOCASSIN ions, it generates datasets for **28 newly available ions in CHIANTI** (e.g., C V, C VI, N VI, N VII, O VII, S I, Cr II, Fe III, Fe IV, Fe IX, Ni II).

---

## 2. Reverse-Engineered File Format Specification

MOCASSIN's `makeElements` subroutine ([`source/hydro_mod.f90`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/source/hydro_mod.f90#L521-L663)) parses each atomic data file using standard list-directed and formatted I/O. Each file consists of **8 sequential sections**:

```
+-------------------------------------------------------------+
| 1. Comments Header                                          |
|    Line 1: NCOMS (Integer: total comment lines)             |
|    Lines 2 .. NCOMS+1: Citations, version, author notes     |
+-------------------------------------------------------------+
| 2. Header Line                                              |
|    NLEVS   NTEMPS (e.g., "10  23")                          |
+-------------------------------------------------------------+
| 3. Level Labels                                             |
|    NLEVS lines, Fortran format (A20)                        |
|    e.g., "   1  2s2 2p 2P0     "                            |
+-------------------------------------------------------------+
| 4. Temperature Grid                                         |
|    NTEMPS lines, format (F10.1), in Kelvin                  |
|    Default: 23 points from 100.0 to 2511886.4 K             |
+-------------------------------------------------------------+
| 5. Mode Flag (iRats)                                        |
|    0 = Collision strengths Upsilon(T) (standard for CHIANTI)|
|    1 = Collision rates / 10**iRats (used only by Fe II)     |
+-------------------------------------------------------------+
| 6. Collision Strength Matrix                                |
|    For each transition (i < j):                             |
|       i   j   Upsilon(T_0)                                  |
|       0   0   Upsilon(T_1)                                  |
|       ...                                                   |
|       0   0   Upsilon(T_{NTEMPS-1})                         |
|    Terminator line: " 0 0 0"                                |
+-------------------------------------------------------------+
| 7. Radiative Transition Probabilities (Einstein A-values)   |
|    Ordered as:                                              |
|       do k = 1, NLEVS - 1                                   |
|          do l = k + 1, NLEVS                                |
|             k   l   A_{lk}  (s^-1)                          |
|    Exactly NLEVS * (NLEVS - 1) / 2 lines                    |
+-------------------------------------------------------------+
| 8. Energy Levels & Statistical Weights                      |
|    do i = 1, NLEVS                                          |
|       i   g_i   E_i (cm^-1, ground level = 0.0)             |
+-------------------------------------------------------------+
```

### Critical Implementation Rules

#### 1. The `if (qx == 0.d0) exit` Rule
In [`source/hydro_mod.f90:609`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/source/hydro_mod.f90#L609):
```fortran
read(121, *) indexlow(2), indexup(2), qx
if (qx == 0.d0) exit
```
Fortran uses `qx == 0.d0` as the loop termination signal. If any collision strength in the table is formatted as `0.00e+00` (which can happen at very high temperatures where spline fits drop to zero), **Fortran exits the collision loop prematurely** and attempts to parse the remaining collision lines as A-values, causing immediate data corruption.

> [!IMPORTANT]
> The generator guarantees that all valid collision strengths are strictly positive by enforcing a floor of `1.0e-30`. Only the dedicated trailer `0 0 0` emits `qx = 0.0`.

#### 2. Array Dimension Limit (`nForLevels = 17`)
In [`source/constants_mod.f90:65`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/source/constants_mod.f90#L65):
```fortran
integer, parameter :: nForLevels = 17           ! number of levels
integer, parameter :: nForLevelsLarge = 142     ! number of levels for FeII
```
And in [`source/emission_mod.f90:1351`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/source/emission_mod.f90#L1351):
```fortran
if (atomic_data_array(elem,ion)%nlevs > size(fLineEm(1,:))) then
   print*, '! equilibrium: model ion has more levels than allowed by nForLevels - please enlarge'
   stop
```
MOCASSIN allocates internal cooling arrays for up to 17 levels. The generator caps all general ions at `NLEVS <= 17` by default.

#### 3. The 259-Entry Mapping in `data/fileNames.dat`
`data/fileNames.dat` lists filenames sequentially corresponding to:
```fortran
do elem = 3, nElements             ! Elements 3 (Li) to 30 (Zn)
   do ion = 1, min(elem+1, 10)     ! Ions I to min(Z+1, X)
      read(17, '(A20)') dataFile(elem, ion)
```
Total entries: $4 + 5 + 6 + 7 + 8 + 9 + 10 + (21 \times 10) = 259$ files. If a file is absent on disk, MOCASSIN marks `lgDataAvailable(elem, ion) = .false.` and skips it.

#### 4. The `feii.dat` Exception (`iRats = 1`)
`feii.dat` uses a dedicated 142-level model from Nahar & Pradhan (1994, 1995) read by a custom branch in `hydro_mod.f90:588` expecting 2,502 blocks of 4 A-values per line. The generator automatically preserves `feii.dat` unless explicitly commanded to overwrite it.

---

## 3. CHIANTI Physics & Numerical Methods

The converter reads raw CHIANTI tables and calculates thermally-averaged collision strengths $\Upsilon(T)$:

### Input File Formats
1. **Energy Levels (`.elvlc`)**: Contains level index, electron configuration, spectroscopic terms, $J$, statistical weight $g = 2J+1$, and energy levels in $\text{cm}^{-1}$ (prefers $E_{\rm obs}$, falling back to $E_{\rm th}$).
2. **Radiative Transitions (`.wgfa`)**: Contains lower level $i$, upper level $j$, transition wavelength, and Einstein A-coefficient $A_{ji}$ in $\text{s}^{-1}$.
3. **Collision Strengths (`.splups` and `.scups`)**:
   - `.splups`: Burgess & Tully (1992) 5-point, 9-point, or multi-point scaled spline fit coefficients.
   - `.scups`: CHIANTI v8+ format providing explicit scaled temperatures `btemp` and scaled collision strengths `bscups`.

### Burgess & Tully (1992) Descaling Formulas
Given transition energy $\Delta E = E_j - E_i$ (Rydberg) and scaling parameter $C$:
$$k_B T / \Delta E = \frac{k_B T}{\Delta E \times 2.17987197 \times 10^{-11}\text{ erg}}$$

The dimensionless scaled temperature $x$ and descaling relation depend on transition type:
- **Type 1 (Electric Dipole Allowed)**:
  $$x = 1 - \frac{\ln C}{\ln(k_B T/\Delta E + C)}, \quad \Upsilon(T) = y(x) \times \ln(k_B T/\Delta E + e)$$
- **Type 2 (Born / Dipole Forbidden)**:
  $$x = \frac{k_B T/\Delta E}{k_B T/\Delta E + C}, \quad \Upsilon(T) = y(x)$$
- **Type 3 (Exchange)**:
  $$x = \frac{k_B T/\Delta E}{k_B T/\Delta E + C}, \quad \Upsilon(T) = \frac{y(x)}{k_B T/\Delta E + 1}$$
- **Type 4 (Resonance / Forbidden)**:
  $$x = 1 - \frac{\ln C}{\ln(k_B T/\Delta E + C)}, \quad \Upsilon(T) = y(x) \times \ln(k_B T/\Delta E + C)$$
- **Type 5 (Dielectronic)**:
  $$x = \frac{k_B T/\Delta E}{k_B T/\Delta E + C}, \quad \Upsilon(T) = \frac{y(x)}{k_B T/\Delta E}$$
- **Type 6 (Proton)**:
  $$x = \frac{k_B T/\Delta E}{k_B T/\Delta E + C}, \quad \Upsilon(T) = y(x)$$

Interpolation $y(x)$ is evaluated via cubic spline (`scipy.interpolate.splrep`/`splev`) on the standard MOCASSIN 23-point logarithmic temperature grid ($T = 10^{2.0} \to 10^{6.4}\text{ K}$, step $0.2\text{ dex}$).

---

## 4. Default Level Counts (`NLEVS`)

To ensure physical consistency with established MOCASSIN nebular cooling calculations while respecting `nForLevels = 17`, the generator includes a built-in mapping of level counts:

| Ion / File | NLEVS | Configurations / Terms Included |
|:---|:---:|:---|
| **C II** (`cii.dat`) | **10** | $2s^2 2p\ ^2P^o$, $2s 2p^2\ ^4P$, $^2D$, $^2S$, $^2P$ |
| **C III** (`ciii.dat`) | **10** | $2s^2\ ^1S$, $2s 2p\ ^3P^o$, $^1P^o$, $2p^2\ ^3P$, $^1D$, $^1S$ |
| **C IV** (`civ.dat`) | **15** | $2s\ ^2S$, $2p\ ^2P^o$, $3s\ ^2S$, $3p\ ^2P^o$, $3d\ ^2D$ |
| **N II** (`nii.dat`) | **15** | Ground configuration $2s^2 2p^2$ and $2s 2p^3$ terms |
| **N III** (`niii.dat`) | **10** | $2s^2 2p\ ^2P^o$, $2s 2p^2\ ^4P$, $^2D$, $^2S$, $^2P$ |
| **O I** (`oi.dat`) | **7** | Complete low-lying atom (fine-structure ground + excited) |
| **O II** (`oii.dat`) | **13** | Ground configuration $2p^3\ ^4S^o, ^2D^o, ^2P^o$ + $2s 2p^4$ |
| **O III** (`oiii.dat`) | **15** | Ground configuration $2p^2\ ^3P, ^1D, ^1S$ + $2s 2p^3$ |
| **Ne II** (`neii.dat`) | **2** | Ground fine-structure doublet $^2P^o_{3/2, 1/2}$ ($12.81\,\mu\text{m}$) |
| **Ne III** (`neiii.dat`)| **9** | $2p^4\ ^3P, ^1D, ^1S$ + $2s 2p^5$ |
| **S II** (`sii.dat`) | **5** | Classic 5-level atom ($^4S^o, ^2D^o_{3/2, 5/2}, ^2P^o_{1/2, 3/2}$) |
| **S III** (`siii.dat`) | **5** | Classic 5-level atom ($^3P_{0, 1, 2}, ^1D_2, ^1S_0$) |
| **Fe VI** (`fevi.dat`) | **17** | Ground term fine structure and low-lying metastable states |
| *Other Ions* | **$\le 17$** | Capped at $\min(\text{total levels}, 15)$ |

---

## 5. Prerequisites & Setup

The tool requires Python 3.9+ with `numpy` and `scipy`:

```bash
# Using the active Miniforge environment:
/home/hmonteiro/miniforge3/bin/python -m pip install numpy scipy
```

---

## 6. CLI Usage & Examples

### Command-Line Arguments Reference

| Option | Type | Default | Description |
|:---|:---:|:---:|:---|
| `--chianti-dir` | Path | None | Path to uncompressed CHIANTI directory (auto-detects `$XUVTOP` or local Cloudy database) |
| `--chianti-tar` | Path | None | Path to a `CHIANTI_*.tar.gz` archive (reads directly without disk extraction) |
| `--download-version` | String | None | Version number (e.g. `10.1`, `11.0.2`) to download from `chiantidatabase.org` |
| `--output-dir` | Path | `data` | Directory where generated `.dat` files will be written |
| `--file-list` | Path | `data/fileNames.dat` | Master ion list mapping |
| `--ions` | String | None | Comma-separated list of ions (e.g. `c_2,o_3,fe_6` or `all`) |
| `--max-levels` | Integer | `17` | Global cap on level count to protect MOCASSIN's array limits |
| `--nlevs` | String | None | Custom per-ion level overrides (e.g. `c_2=10,o_3=15`) |
| `--version-tag` | String | Auto | Version string in comments header |
| `--skip-existing` | Flag | False | Skip files that already exist in output directory |
| `--compare-with` | Path | None | Directory with reference files to print line-count diffs |

---

### Common Workflows

#### Workflow 1: Generate All Available Ions from Local CHIANTI Data
```bash
/home/hmonteiro/miniforge3/bin/python tools/chianti_converter/generate_chianti_atomic_data.py \
    --chianti-dir /home/hmonteiro/software/cloudy-c25.00/data/chianti \
    --output-dir data.chianti-generated
```

#### Workflow 2: Read Directly from a Downloaded CHIANTI Archive
```bash
/home/hmonteiro/miniforge3/bin/python tools/chianti_converter/generate_chianti_atomic_data.py \
    --chianti-tar ~/Downloads/CHIANTI_10.1_database.tar.gz \
    --output-dir data
```

#### Workflow 3: Automatically Download Official CHIANTI 10.1 Release
```bash
/home/hmonteiro/miniforge3/bin/python tools/chianti_converter/generate_chianti_atomic_data.py \
    --download-version 10.1 \
    --output-dir data.chianti-10.1
```

#### Workflow 4: Update Specific Ions with Custom Level Counts
```bash
/home/hmonteiro/miniforge3/bin/python tools/chianti_converter/generate_chianti_atomic_data.py \
    --ions c_2,c_3,o_3,s_2 \
    --nlevs c_2=10,o_3=15,s_2=5 \
    --output-dir data \
    --compare-with data.chianty-10
```

---

## 7. Verification & Automated Test Suite

A dedicated regression test suite is provided in [`test_converter.py`](file:///home/hmonteiro/software/mocassin_HM_2026_3DPDR_testing/tools/chianti_converter/test_converter.py):

```bash
/home/hmonteiro/miniforge3/bin/python tools/chianti_converter/test_converter.py
```

### Verification Matrix
- **Test 1: Burgess–Tully Descaling**: Checks mathematical exactness for Type 1, 2, 3, and 4 transitions against analytical bounds.
- **Test 2: C II Benchmark Comparison**: Compares newly generated `cii.dat` against `data.chianty-10/cii.dat`:
  - 45/45 Einstein A-coefficients match to full precision.
  - 10/10 energy levels and statistical weights match.
  - Median relative difference in $\Upsilon(T)$ across all 1,035 matrix values over 5 orders of magnitude in temperature is **$< 0.76\%$**.
- **Test 3: Fortran Reader Runtime Compatibility**: Compiles an exact standalone `gfortran` test harness mirroring `source/hydro_mod.f90` and validates that **all 124 generated atomic data files parse with 0 exit code and 0 runtime errors**.
- **Test 4: Premature Termination Guard**: Verifies that collision strength blocks never contain premature `0.00e+00` values that would trigger early Fortran loop exit.
- **Test 5: Official CHIANTI Directory Resolution & .scups Support**: Verifies automated case-insensitive directory resolution (e.g. `CHIANTI_10.1` -> `chianti_v10.1`) and parses standard CHIANTI 10.1 `.scups` files.
