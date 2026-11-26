# Deep Forensic Audit: `lgCylinder2D` & Planar 2D Geometry, Volumes, Averaging, and Packet Transfer

- **Repository:** `/home/hmonteiro/software/mocassin_HM_2026`
- **Branch:** `feature/2D-mode`
- **Baseline Audited Commit:** `67202ee` (*Fix 2D writeTauNu O(N_max) hang and eliminate duplicate MPI execution*)
- **Audit Date:** 2026-11-25
- **Scope:** Full source audit across [`source/*.f90`](file:///home/hmonteiro/software/mocassin_HM_2026/source) covering `lgCylinder2D` (`cylinder`) and planar 2D closures (`noSpecular`, `specular`, `buffer`) across three pillars:
  1. **Volumes** (`getVolume`, `getVolumeCon`, mass integrals, serialization, light echoes)
  2. **Averaging & Normalization** (`JSte`/`JDif` MC estimators, `outputGas`, `Tdust`, `dustHeatingBudget`, `slit`, `writeSED`)
  3. **Packet Transfer** (`photon_mod`, `emission_mod`, `pathIntegration_mod`, `bufferTransfer2D`, `fluorescence_mod`, subgrids, sources)

> **Executive Verdict:** For standard **axisymmetric equatorial-symmetric models (`2D XZ cylinder` + `symmetricXYZ`, single on-axis star, no `resLinesTransfer` ray-tracing, no warm restart)**, the cylindrical ring volumes, Monte Carlo $J_\nu$ estimator normalization, dust radiative equilibrium temperatures, and volume-integrated line/ionic averages are **mathematically exact**. Outside that baseline path, this audit identifies **16 tracked issues** (3 Critical, 4 High, 5 Medium, 4 Low/Unsupported) listed below so fixes can be implemented and checked off systematically.

---

## Master Fix-Tracking Checklist

Use this table to track remediation commits on `feature/2D-mode` (change `[ ]` to `[x]` and record the fixing commit hash).

| Status | ID | Pillar | Severity | Primary Files & Lines | Summary | Fix Commit / Notes |
| :---: | :---: | :--- | :--- | :--- | :--- | :--- |
| `[x]` | **P0** | Transfer | Critical | `iteration_mod.f90:949`, `output_mod.f90:2595`, `pathIntegration_mod.f90:15-363` | `writeTauNu` executed by all MPI ranks & `(0,1,0)` perpendicular bounce loop (`maxTau=1e7`) in planar/cylinder 2D | Fixed in **`67202ee`** (verified ~360-460x speedup) |
| `[x]` | **T1** | Transfer | **Critical** | `emission_mod.f90:1940-2542` | `setResLineEscapeProb` has no `lg2D` / `lgCylinder2D` support; rays abort on step 1 across `ny=1`, breaking dusty Ly$\alpha$/resonance heating | Fixed in **`187c168`**: added cylindrical & planar 2D ray tracing (dusty PN150 cyl vs 3D: $H\beta$ $-1.44\%$, $T_e$ $+0.08\%$, Ly$\alpha$ $0.5139$ vs $0.5433$) |
| `[x]` | **A2** | Averaging | **Critical** | `iteration_mod.f90:93-110,692-706`, `hydro_mod.f90:674-688`, `emission_mod.f90:833-851`, `photon_mod.f90:715-759`, `grid_mod.f90` | Out-of-bounds write on `JDif(:, 1:nbins)` inside `do freq = 1, nLines` when `lg2D .and. lgDebug`, dead `TwoDscaleJ`, and upstream `linePDF`/`nLines`/`getNu2` bounds mismatches | Fixed in **`7818329`**: removed `TwoDscaleJ` & invalid `JDif` loop; lazy-allocated `linePDF`/`linePackets`; aligned `nLines`/`linePDF` (`175` HI, `34` HeI) and `getNu2` discrete CDF bounds (`-fbounds-check`: $H\beta$ $-0.01\%$ vs no-debug, `[OIII]` An/MC $=1.0125$; gate PASS) |
| `[x]` | **V6** | Volumes / Restart | **Critical** | `grid_mod.f90:3158-3162,3369-3500,3647-3685`, `output_mod.f90:440`, `mocassinPlot.f90:172,213`, `mocassinWarm.f90:55-58` | `grid3.out` omits `lgCylinder2D`, `lgBuffer2D`, `lgSpecular2D`; `mocassinWarm`/`Output`/`Plot` reset `lgCylinder2D = .false.` (+ pre-read loop in `mocassinWarm`) | Fixed in **`7c6e580`**: serialized `lgCylinder2D/lgBuffer2D/lgSpecular2D` in `grid3.out` with legacy `iostat` fallback; extracted `getWeight2D`; fixed `mocassinWarm` pre-read & `mocassinPlot` `dV` (`mocassinOutput` $0.000\%$ diff, `grid4.out` $0.0000\%$ $V_{\mathrm{oct}}$ err, `mocassinWarm` $H\beta$ $-0.10\%$; gate PASS) |
| `[x]` | **V2** | Volumes / Transfer | **High** | `grid_mod.f90:600-625,1528-1615`, `pathIntegration_mod.f90:36` | `lgCylinder2D` without `symmetricXYZ` puts `xAxis` on `[-Rnx,+Rnx]` instead of `[0,+Rnx]`, yielding **negative ring volumes** for $x < 0$ | Fixed in **`ce63993`**: automatic `lgCylinder2D` without `symmetricXYZ` maps `xAxis` on `[0,Rnx]` & `zAxis` on `[-Rnz,+Rnz]` (odd `nz`); 3D `[-R,+R]^3` `densityFile` slices $R \ge 0$ half-plane (`zP==0` bottom clamp fixed in `initPathTauStart`); auto nosym $H\beta$ $-0.00\%$ & `dfile` nosym $H\beta$ $-0.05\%$; gate PASS |
| `[x]` | **T2** | Transfer | **High** | `photon_mod.f90:1462-1511,1727-1740` | Grazing/on-wall `abs(dSx) < 1e-10` or `abs(dSz) < 1e-10` replaces `0` with a full mean grid cell width instead of advancing `xP`/`zP` | Fixed in **`b72a89a`**: clamped `dSz >= 0`, allowed zero-distance `C_out >= 0` / `C_in <= 0` wall advance (`B_cyl` tangency guard `-1e-6 * sqrt(R2*A)`), removed mean-cell-width substitution & `1e-10` floor; stretched `(R,Z)` HII40 $H\beta$ $+0.12\%$, $T_e$ $+0.09\%$, `[OIII]` $+0.76\%$ vs 3D; gate PASS |
| `[x]` | **T3** | Transfer | **High** | `photon_mod.f90:2348-2525,3715-3915` | Planar `bufferTransfer2D` ignores `totalLines` (line escape cooling) and all dust scattering/`dustPDF` (crashes if `noGas` + `buffer`) | Fixed in **`db8f730`**: implemented dust scattering (`scaOpac/opacity` with HG/isotropic), pure-dust `dustPDF` thermal emission (`noGas`), and gas `totalLines` optically thin line escape (`if (enPacket%lgLine) return` at all 6 call sites); verified `noGas` (`tau1.000` + `p0tau1` sphere `buffer`) + `HII40`/`PN150` `buffer` ($\le 0.26\%$); gate PASS |
| `[x]` | **V1** | Volumes | **High** | `continuum_mod.f90:573-636` | `getVolumeCon` lacks `lgCylinder2D` branch, corrupting `setLdiffuse` weights `LdiffuseLoc` by $\propto \delta_{\rm perp}/R$ | Fixed in **`5ed4207`**: added exact `lgCylinder2D` toroidal volume branch (`R_in_c`, `R_out_c`, `dz`, octant/full factor) in `getVolumeCon`; verified `diffuseSource 15.0` on PN75 (`2D_cyl` vs `3D`: $T_e$ $+0.17\%$, `[OIII]` $+0.37\%$, $H\beta$ $-3.62\%$ discrete coarse sphere); gate PASS |
| `[x]` | **V4** | Volumes | **Medium** | `grid_mod.f90:3342-3412,3482-3540`, `continuum_mod.f90:573,635-730` | Planar 2D outer-edge cells (`xP==nx`, `yP==ny`, `zP==nz`) omit `/ 2.0` when `lgSymmetricXYZ=.true.`, making boundary volumes $2\times$ ($4\times$ at corners) too large | Fixed in **`be486c5`**: added `/ factor` (`2.0` if `lgSymmetricXYZ`) at outer boundaries in planar 2D `getVolume` and planar 2D/3D `getVolumeCon`, plus outer half-cell centroid in `getWeight2D`; verified exact octant box cylinder volume (`0.000011%` error in `grid4.out`), `HII40_2D_buffer` (`Hb=14.8914`), and 4-model gate ($\le 0.64\%$) |
| `[x]` | **V5** | Volumes | **Medium** | `grid_mod.f90:755,811,1991,2868` | Planar 2D `totalGasMass`, `totalDustMass`, and `inputGasMass`/`inputDustMass` scaling (`massFac`) omit revolution weight `weight2D` | Fixed in **`f04a3d1`**: weighted `dV = getVolume(...) * getWeight2D(...)` in `fillGrid` (`gasMass`/`dustMass` loops), `setMotherGrid`, and `setSubGrids`; verified `HII40` planar `totalGasMass` & `gasMass` scaled `Hden` match `2D_cylinder` to `0.00%` (`+3.50%` vs `3D`); gate PASS ($\le 1.29\%$) |
| `[x]` | **V3** | Volumes / Setup | **Medium** | `set_input_mod.f90:888-922`, `grid_mod.f90:3311`, `continuum_mod.f90:606` | `readInput` allows `2D XY cylinder` / `2D YZ cylinder`, but `getVolume`/`photon_mod`/`pathIntegration_mod` hardcode `xAxis`=$R$, `zAxis`=$Z$, `yP=1` | **Fixed**: enforced `(XZ, Z)` orientation and `nx >= 2, nz >= 2` for `lgCylinder2D` in `set_input_mod.f90` (and closed the `else if (lg2D)` bypass so `planeIonization` and downstream checks execute in 2D), plus `nx, nz < 2` runtime guards in `getVolume`/`getVolumeCon`; verified clean rejection of `2D XY cylinder`, `nx=1` cylinder, and `planeIonization+2D`, valid `2D XY Y buffer` PASS; gate PASS ($\le 1.00\%$) |
| `[ ]` | **A1** | Averaging | **Medium** | `update_mod.f90:2497-2505`, `iteration_mod.f90:1032-1047` | `dustHeatingBudget` sums over cells without multiplying by `dV`, biasing dust heating budget toward $R \to 0$ | _Open_ |
| `[ ]` | **A3** | Averaging | **Medium** | `output_mod.f90:318-333`, `mocassinPlot.f90:177-190` | `slit` filter checks Cartesian point `yAxis(1) = 0` instead of arc/chord fraction of torus $R_i$ inside slit aperture | _Open_ |
| `[ ]` | **A4** | Averaging | **Low** | `output_mod.f90:2986-3092` | `writeSED` general/2D caveats: `ijk=0` cavity continuum skipped, gas-only `SEDratio` zero-division, planar 2D perp-wall SED distortion | _Open_ |
| `[ ]` | **T4** | Transfer | **Low** | `photon_mod.f90:161-163,224-226,878-880` | `"diffExt"` `posDiff%z` uses `cellLoc(2)` instead of `cellLoc(3)` + discrete `(xAxis(ix), 0, zAxis(iz))` launch for `diffExt`/`resLinePackets` | _Open_ |
| `[ ]` | **V7 / T5 / T6** | Guards | **Low** | `grid_mod.f90:4056-4260`, `photon_mod.f90:2166-2283`, `fluorescence_mod.f90:423-860` | Missing `lgCylinder2D` support or startup guards for `vEcho`, `nGrids > 1`, off-axis stars, and `fluorescence_mod` | _Open_ |

---

## 1. Volumes Audit

### 1.1 Verified Behavior (`lgCylinder2D` + `symmetricXYZ` in `XZ`)

1. **Toroidal Cell Volume (`getVolume` in [grid_mod.f90:3209-3240](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3209-L3240)):**
   - Radial walls: $R_{\mathrm{in},1} = 0$, $R_{\mathrm{out},i} = (x_{i+1} + x_i)/2$ for $1 \le i < n_x$, and $R_{\mathrm{out},n_x} = x_{n_x}$ (matching the photon outer boundary in `photon_mod.f90:1414-1423`).
   - Axial heights: $\Delta Z_1 = (z_2 - z_1)/2$, $\Delta Z_k = (z_{k+1} - z_{k-1})/2$ for $1 < k < n_z$, and $\Delta Z_{n_z} = (z_{n_z} - z_{n_z-1})/2$ (matching `Z_in`/`Z_out_c` in `photon_mod.f90:1426-1435`).
   - Volume in units of $10^{45}\,\mathrm{cm^3}$ when `lgSymmetricXYZ = .true.`:
     $$V_{\mathrm{oct}}(i, k) = \frac{\pi}{4}\left(R_{\mathrm{out},i}^2 - R_{\mathrm{in},i}^2\right)\Delta Z_k \times 10^{-45}$$
2. **Full-Nebula Mass & Emission Integrals:**
   - Multiplying by `8.` in `fillGrid` ([grid_mod.f90:851-860](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L851-L860)) and `outputGas` ([output_mod.f90:1041-1043](file:///home/hmonteiro/software/mocassin_HM_2026/source/output_mod.f90#L1041-L1043)) yields
     $$8\,V_{\mathrm{oct}}(i, k) = 2\pi\left(R_{\mathrm{out},i}^2 - R_{\mathrm{in},i}^2\right)\Delta Z_k = \int_0^{2\pi} d\phi \int_{R_{\mathrm{in},i}}^{R_{\mathrm{out},i}} R\,dR \,(2\,\Delta Z_k)$$
     covering the full $2\pi$ azimuth across both $+Z$ and $-Z$ hemispheres. Thus `totalGasMass`, `totalDustMass`, `totalSpecMass`, `inputGasMass`/`inputDustMass` density scaling (`massFac`), and `HbetaVol` are **exact** in symmetric `2D XZ cylinder` runs.
3. **Volume Cache (`initCellVolumes` in [photon_mod.f90:3902-3926](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L3902-L3926)):**
   - Populates `cell_volumes(ig)%v(ix,iy,iz) = getVolume(grid(ig), ix, iy, iz)` accurately for all active cells.

---

### 1.2 Detailed Volume Findings

#### `[x]` **V1 (High) — Missing `lgCylinder2D` Branch in `getVolumeCon` / `setLdiffuse`**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): added the exact `lgCylinder2D` toroidal ring volume calculation inside `getVolumeCon` (`continuum_mod.f90:573-636`, matching `getVolume` in `grid_mod.f90`). Verified on `PN75` with dominant `diffuseSource 15.0 75000. "blackbody"` in `2D XZ Z cylinder` vs `3D`: `<Te(H+)>` $+0.17\%$, `[OIII] 5008` $+0.37\%$, `HeII 4686` $+3.25\%$, `Hbeta` $-3.62\%$ (consistent with $13^3$ vs $13\times 13$ discrete spherical boundary), and 4-model regression gate ($\le 0.96\%$).
- **Location:** [`source/continuum_mod.f90:562-774`](file:///home/hmonteiro/software/mocassin_HM_2026/source/continuum_mod.f90#L562-L774) (`getVolumeCon`), called at [`source/continuum_mod.f90:498,514,531,546`](file:///home/hmonteiro/software/mocassin_HM_2026/source/continuum_mod.f90#L498-L546) (`setLdiffuse`).
- **Root Cause:** Under `else if (lg2D)` (`line 602`), `getVolumeCon` only computes the flat Cartesian slab `dx * dy * dz` with `dy = deltaPerp2D`, without checking `lgCylinder2D`.
- **Impact:** When a `2D cylinder` model specifies `Ldiffuse > 0` (diffuse extra continuum source), `setLdiffuse` distributes `LdiffuseLoc(iloc) = Ldiffuse * (Hden * dV) / sum(Hden * dV)` using flat slab `dV` instead of cylindrical ring volume $\propto (R_{\mathrm{out}}^2 - R_{\mathrm{in}}^2)\Delta Z$. Inner axial cells ($R \to 0$) receive an artificially huge share of `Ldiffuse` relative to outer rings by factor $\propto \delta_{\mathrm{perp}} / R$.
- **Recommended Fix:** Make `getVolumeCon` call or match `getVolume` (e.g., `use grid_mod, only : getVolume` is circular if `grid_mod` uses `continuum_mod` — check module dependency order; if circular, replicate the exact `lgCylinder2D` block and `/ 2.0` outer-wall fix in `getVolumeCon`).

#### `[x]` **V2 (High) — `lgCylinder2D` Without `symmetricXYZ` Generates Negative Ring Volumes**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated `fillGrid` (`case ("Y")`) when `lgCylinder2D` without `lgSymmetricXYZ` to generate `xAxis` over `[0, Rnx]` and `zAxis` over `[-Rnz, +Rnz]` (odd `nz`); added $R \ge 0$ half-plane slicing in `setMotherGrid` for 3D `[-R, +R]^3` `densityFile` inputs; fixed `zP == 0` at `zAxis(1)` in `initPathTauStart` (`pathIntegration_mod.f90`). Verified on `PN75` `2D XZ Z cylinder` without `symmetricXYZ` (`Hb=0.05396`, `-0.00%` vs symmetric) and 3D `25^3` `densityFile` (`Hb=0.05393`, `-0.05%`), plus 4-model regression gate ($\le 1.15\%$).
- **Location:** [`source/grid_mod.f90:600-625,1528-1615`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L600-L625), [`source/pathIntegration_mod.f90:36`](file:///home/hmonteiro/software/mocassin_HM_2026/source/pathIntegration_mod.f90#L36).
- **Root Cause:**
  1. In `getVolume` (`line 3238-3239`), when `.not. lgSymmetricXYZ`, `lgCylinder2D` computes `Pi * (R_out^2 - R_in^2) * dz` (full $2\pi$ azimuthal ring over $[-Z_{\max}, +Z_{\max}]$), which is algebraically valid **if `xAxis` spans $R \in [0, R_{\max}]$ and `zAxis` spans $Z \in [-Z_{\max}, +Z_{\max}]$**.
  2. However, `fillGrid` (`line 606`) initializes `xAxis` over `[-Rnx, +Rnx]`, and `setMotherGrid` (`line 1224`) reads all `x3D` values (`[-Rmax, +Rmax]`) from a 3D density file!
  3. Furthermore, at `zP == 1` (`line 3223`) and `zP == nz` (`line 3227`), `getVolume` sets `dz = (zAxis(2)-zAxis(1))/2.` matching photon exit at `zAxis(1)` (`line 1427`, `2977`), **but** for `xAxis` starting at `-Rnx`: $R_{\mathrm{out},i}^2 - R_{\mathrm{in},i}^2 < 0$ for all $x < 0$, and at `i = 1`, `getVolume` hardcodes `R_in = 0.0`, giving `Pi * ((-Rnx + dx/2)^2 - 0) * dz > 0` instead of a thin outer ring!
- **Impact:** Any `2D cylinder` simulation with `.not. lgSymmetricXYZ` (north–south asymmetric bipolar/disk models) produces negative volumes and broken $R_{\mathrm{cyl}} \ge 0$ indexing.
- **Recommended Fix:** In `lgCylinder2D`, `xAxis` represents cylindrical radius $R \ge 0$ **regardless of `lgSymmetricXYZ`**:
  - In `fillGrid`, always populate `xAxis(1:nx)` in `[0, Rnx]` when `lgCylinder2D` is `.true.`, while `zAxis(1:nz)` spans `[0, Rnz]` if `lgSymmetricXYZ` or `[-Rnz, +Rnz]` (odd `nz`) if `.not. lgSymmetricXYZ`.
  - When `lgDfile` + `lgCylinder2D` + `.not. lgSymmetricXYZ`, validate or extract only the $R = x \ge 0$ half-plane (or adjust `nx` accordingly), and enforce `xAxis(1) >= 0.0` at startup. Also handle nonzero inner-radius grids where `xAxis(1) > 0` by checking whether `R_in` at `xP == 1` should be `0` or `xAxis(1)` (if `xAxis(1) > 0`, `R_in = xAxis(1)` matches packet boundaries? Note `R_in = 0` at `xP = 1` assumes `xAxis(1) = 0`, which holds in auto grids).

#### `[x]` **V3 (Medium) — Hardcoded `(XZ, Z)` Orientation for `lgCylinder2D`**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated `readInput` (`set_input_mod.f90:888-922`) to enforce `plane2D == "XZ"` and `symmetryAxis2D == "Z"` plus `nxIn(1) >= 2 .and. nzIn(1) >= 2` when `lgCylinder2D`, closed the `else if (lg2D)` validation branch before `if (lgPlaneIonization .and. lg2d)` so downstream compatibility checks execute in 2D mode, and added `nx < 2 .or. nz < 2` runtime guards in `getVolume` (`grid_mod.f90:3311`) and `getVolumeCon` (`continuum_mod.f90:606`). Verified clean diagnostic termination on `2D XY Y cylinder`, `nx=1` cylinder, and `planeIonization` + 2D, successful execution of valid planar `2D XY Y buffer`, and 4-model regression gate ($\le 1.00\%$).
- **Location:** [`source/set_input_mod.f90:212-308`](file:///home/hmonteiro/software/mocassin_HM_2026/source/set_input_mod.f90#L212-L308), [`source/grid_mod.f90:3209-3240`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3209-L3240), [`source/photon_mod.f90:1409-1486`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L1409-L1486).
- **Root Cause:** `readInput` allows `2D XY ... cylinder` (`nz = 1`, `perpAxis2D = "Z"`) and `2D YZ ... cylinder` (`nx = 1`, `perpAxis2D = "X"`), yet `getVolume`, `photon_mod`, and `pathIntegration_mod` assume `perpAxis2D == "Y"`, `ny = 1`, `xAxis = R`, and `zAxis = Z` (`nx > 1`, `nz > 1`). Choosing `2D XY cylinder` accesses `grid%zAxis(2)` out-of-bounds.
- **Recommended Fix:** In `set_input_mod.f90`, reject `lgCylinder2D` unless `plane2D == "XZ"` (and `symmetryAxis2D == "Z"`), printing an informative error and stopping (or dynamically map `(R, Z)` axes to the chosen plane). Also guard `grid%nx < 2` or `grid%nz < 2` before indexing `xAxis(2)` / `zAxis(2)` in `getVolume`.

#### `[x]` **V4 (Medium) — Planar 2D Outer Boundary Volumes $2\times$ ($4\times$) Too Large Under `lgSymmetricXYZ`**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated planar 2D `getVolume` (`grid_mod.f90:3342-3412`) and planar 2D/3D `getVolumeCon` (`continuum_mod.f90:573,635-730`) so outer boundary cells (`xP == nx`, `yP == ny`, `zP == nz`) divide by `factor` (`2.0` when `lgSymmetricXYZ`, `1.0` otherwise), and updated `getWeight2D` (`grid_mod.f90:3482-3540`) to evaluate outer boundary rings at their half-cell centroid `axis(n) - abs(axis(n) - axis(n-1)) / 4.0` under `lgSymmetricXYZ`. Verified via `mocassinPlot` `grid4.out` unclipped `13x1x13` planar `2D XZ Z buffer` octant volume sum (`2.155133e12` vs exact $(\pi/4) R^2 Z = 2.155133\times 10^{12}$, error `0.000011%`), `HII40_2D_buffer` ($H\beta = 14.8914$), and 4-model regression gate ($\le 0.64\%$).
- **Location:** [`source/grid_mod.f90:3242-3327`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3242-L3327) (`getVolume` planar branch) and [`source/continuum_mod.f90:602-724`](file:///home/hmonteiro/software/mocassin_HM_2026/source/continuum_mod.f90#L602-L724) (`getVolumeCon`).
- **Root Cause:** Compare planar 2D `getVolume` (`lines 3252-3254, 3264-3266, 3278-3280, 3291-3293, 3306-3308, 3318-3320`) against 3D `getVolume` ([`grid_mod.f90:3331-3362`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3331-L3362)):
  - In 3D, with `lgSymmetricXYZ = .true.`, `factor = 2.0` divides **both** `xP == 1` and `xP == nx` by `2.0`, because photons escape at `xAxis(nx)` (half-cell width `(xAxis(nx) - xAxis(nx-1))/2`).
  - In planar 2D, `xP == 1` divides by `2.0`, **but `xP == nx`, `yP == ny`, `zP == nz` omit `/ 2.0`**, giving full-cell width at the outer wall! (The same omission exists in `getVolumeCon` for both 2D and 3D.)
- **Impact:** Outer boundary cells in planar 2D modes (`noSpecular`, `specular`, `buffer`) have $2\times$ inflated volumes ($4\times$ at outer corners such as `(nx, nz)`), depressing `JSte` by $2\times$–$4\times$ in the outer shell and biasing `outputGas` boundary emission by $>2\times$.
- **Recommended Fix:** Divide outer-edge cell widths `(xAxis(nx)-xAxis(nx-1))`, `(yAxis(ny)-yAxis(ny-1))`, `(zAxis(nz)-zAxis(nz-1))` by `factor` (`2.0` when `lgSymmetricXYZ`, `1.0` otherwise — wait: verify whether 3D `.not. lgSymmetricXYZ` also exits at `xAxis(nx)`! In 3D `.not. lgSymmetricXYZ`, `photon_mod.f90:1535` sets wall at `xAxis(nx)`, so physically outer boundary cells are half-width `(x_nx - x_{nx-1})/2` there too, though keeping 3D compatibility `/ factor` matches existing 3D `getVolume` at `3342-3362`).

#### `[x]` **V5 (Medium) — Planar 2D `totalGasMass`, `totalDustMass`, and `inputGasMass`/`inputDustMass` Scaling Omit `weight2D`**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated `dV = getVolume(...) * getWeight2D(...)` in `fillGrid` (`grid_mod.f90:755,811`, covering `inputGasMass` and `inputDustMass` density scaling), `setMotherGrid` (`grid_mod.f90:1991`), and `setSubGrids` (`grid_mod.f90:2868`). Verified on `HII40` that planar `2D XZ Z buffer` `totalGasMass` (`3.6195e-10`) and `gasMass 100.` scaled `Hden` (`6.8685`) match `2D XZ Z cylinder` to **`0.00%`** (`+3.50%` vs `3D`, matching spherical discretization), and 4-model regression gate passes ($\le 1.29\%$).
- **Location:** [`source/grid_mod.f90:734-868`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L734-L868) (`fillGrid`), [`source/grid_mod.f90:1883-2020`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L1883-L2020) (`setMotherGrid`).
- **Root Cause:** In planar 2D (`lg2D .and. .not. lgCylinder2D`), `dV = getVolume(...)` is only the thin 2D slice thickness `deltaPerp2D * dx * dz`. Neither `setMotherGrid` (`lines 1891-1986`) nor `fillGrid` (`lines 741, 797`) applies the azimuthal revolution factor `weight2D = (Pi/2) * rCyl / deltaPerp2D` (as `outputGas` does at [`output_mod.f90:510-585`](file:///home/hmonteiro/software/mocassin_HM_2026/source/output_mod.f90#L510-L585)).
- **Impact:** Printed total gas/dust masses in planar 2D are only the slice mass rather than the 3D nebula mass, and specifying `gasMass` or `dustMass` in `input.in` scales `Hden`/`Ndust` up by $V_{\mathrm{3D}} / V_{\mathrm{slice}}$ (~$\sim 50\times$). Note: **`lgCylinder2D` is NOT affected** because `getVolume` returns $V_{\mathrm{oct}}$ directly.
- **Recommended Fix:** Apply the shared cylindrical revolution weight `weight2D` (extract a helper in `grid_mod` shared by `fillGrid`, `setMotherGrid`, `outputGas`, and `mocassinPlot`) whenever integrating 3D nebula mass/luminosity in planar 2D mode.

#### `[x]` **V6 (Critical) — `grid3.out` Drops `lgCylinder2D`, `lgBuffer2D`, `lgSpecular2D` (`mocassinWarm`, `mocassinOutput`, `mocassinPlot`)**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): serialized `lgCylinder2D`, `lgBuffer2D`, and `lgSpecular2D` in `writeGrid` (`grid3.out`) and parsed them in `resetGrid` with backward-compatible `iostat` fallback plus `iOrigin/jOrigin/kOrigin = 1` perpendicular pinning; extracted shared `getWeight2D` in `grid_mod.f90` for `output_mod.f90` and `mocassinPlot.f90`; removed spurious pre-read `resetGrid` loop in `mocassinWarm.f90`. Verified on `PN75_2D_cylinder`: `mocassinOutput` matches direct output to `0.000%`, `mocassinPlot` `grid4.out` cylindrical `V_oct` matches theory to `0.0000%`, `mocassinWarm` continues to `0.05390` (`-0.10%`), and 4-model gate PASSES ($\le 0.51\%$).
- **Location:** [`source/grid_mod.f90:3158-3162,3369-3500,3647-3685`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3158-L3162), [`source/mocassinWarm.f90:55-58`](file:///home/hmonteiro/software/mocassin_HM_2026/source/mocassinWarm.f90#L55-L58), [`source/mocassinOutput.f90:44-45`](file:///home/hmonteiro/software/mocassin_HM_2026/source/mocassinOutput.f90#L44-L45), [`source/mocassinPlot.f90:73,172,213`](file:///home/hmonteiro/software/mocassin_HM_2026/source/mocassinPlot.f90#L73).
- **Root Cause:**
  1. `writeGrid` writes only `lg2D, plane2D, symmetryAxis2D, minPerp2D, maxPerp2D, deltaPerp2D` to `output/grid3.out` — omitting `lgCylinder2D`, `lgBuffer2D`, and `lgSpecular2D`.
  2. Standalone executables (`mocassinWarm`, `mocassinOutput`, `mocassinPlot`) call `resetGrid` without calling `readInput()`. Consequently, `lgCylinder2D`, `lgBuffer2D`, and `lgSpecular2D` stay `.false.` (their module defaults in `common_mod`), silently reverting `getVolume` and `photon_mod` to planar `noSpecular`!
  3. Additionally, [`mocassinWarm.f90:56-59`](file:///home/hmonteiro/software/mocassin_HM_2026/source/mocassinWarm.f90#L56-L59) has a pre-existing bug calling `call resetGrid(grid3D(iGrid))` in a loop before `nGrids` is read, followed by `call resetGrid(grid3D)`, and `mocassinPlot.f90:172,213` lacks `weight2D` in planar 2D and projects `2D cylinder` ring luminosities `Hden * dV` only onto `yAxis(1) = 0`.
- **Recommended Fix:** Serialize `lgCylinder2D`, `lgBuffer2D`, `lgSpecular2D` in `writeGrid` (`grid3.out`) and parse them with backward-compatible `iostat` fallback in `resetGrid`; apply `weight2D` in `mocassinPlot` (and clearly document 2D projection meaning for `mocassinPlot` / `writeContCube`); remove/clean the spurious pre-read loop in `mocassinWarm.f90`.

#### `[ ]` **V7 (Low) — Light Echo `vEcho` Incompatible with 2D Grids**
- **Location:** [`source/grid_mod.f90:4175-4260`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L4175-L4260) (`vEcho`), called at [`source/grid_mod.f90:1908-1911`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L1908-L1911).
- **Root Cause:** `vEcho` unconditionally accesses `grid%yAxis(yP+1)` when `yP == 1 == ny`, going out-of-bounds in `2D XZ` (`ny == 1`), and integrates rectangular subcells `(xx, yy)` rather than toroidal light-echo geometry.
- **Recommended Fix:** Reject `lg2D` + `lgEcho` cleanly at input/startup with an explicit diagnostic message (or implement 2D axisymmetric echo integration).

---

## 2. Averaging & Normalization Audit

### 2.1 First-Principles Proof of `JSte`/`JDif` & `outputGas` Normalization in `lgCylinder2D`

1. **Mean Intensity Estimator (`JSte`, `JDif`):**
   - Each packet carries energy $\delta E_i = L_{\star,i} / N_{\mathrm{phot},i}$ ([grid_mod.f90:3434](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L3434)). Thus the launched packet ensemble injects the **full $4\pi$ luminosity $L_\star$**.
   - In `lgCylinder2D` with `lgSymmetricXYZ = .true.`:
     - `initPhotonPacket` ([photon_mod.f90:690-698](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L690-L698)) folds initial stellar directions into the first octant ($+x, +y, +z$).
     - During transfer ([photon_mod.f90:857-860,985-988,2085-2088,3245-3265](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L2085-L2088)), scattered/re-emitted packets roam freely over the **full $2\pi$ azimuth** in $(x, y)$ and reflect **only** at $z = 0$ (equatorial mirror plane), remaining in the upper half-space ($z \ge 0$, combining 4 octants into one ring).
   - Over the 3D sphere, two symmetry-equivalent hemispheres ($\pm Z$, each of $2\pi$ azimuthal ring volume $V_{\mathrm{hemi}} = \pi(R_{\mathrm{out}}^2 - R_{\mathrm{in}}^2)\Delta Z = 4\,V_{\mathrm{oct}}$) combine to form the full two-hemisphere torus $V_{\mathrm{torus}} = 2\pi(R_{\mathrm{out}}^2 - R_{\mathrm{in}}^2)\Delta Z = 8\,V_{\mathrm{oct}}$.
   - Because each photon path length $\ell$ inside $(i, k)$ is summed over all $N_{\mathrm{phot}}$ packets representing total $4\pi$ luminosity $L_\star$, the true $4\pi$-averaged mean intensity in ring $(i, k)$ is:
     $$4\pi J_\nu(i, k) = \frac{\sum \ell\,\delta E}{V_{\mathrm{torus}}(i, k)} = \frac{\sum \ell\,\delta E}{8\,V_{\mathrm{oct}}(i, k)}$$
   - In `pathSegment` ([photon_mod.f90:1864-1875,2146-2157](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L2146-L2157)), the code accumulates $\ell\,\delta E / dV = \ell\,\delta E / V_{\mathrm{oct}}$, and in `iterateMC` ([iteration_mod.f90:712-729](file:///home/hmonteiro/software/mocassin_HM_2026/source/iteration_mod.f90#L712-L729)), `lgSymmetricXYZ = .true.` divides `JSte`, `JDif`, and `escapedPackets` by **`8.`**. The factors cancel **identically**.
   - *(Note on `lgSymmetricXYZ = .false.` once Bug V2 is fixed: `getVolume` returns $V_{\mathrm{hemi}} = \pi(R_{\mathrm{out}}^2 - R_{\mathrm{in}}^2)\Delta Z$ spanning all $\pm Z$, and `iteration_mod.f90:724` does NOT divide by `8.`, so $\sum \ell\,\delta E / V_{\mathrm{hemi}}$ is equally exact.)*
2. **Volume-Weighted Gas Outputs (`outputGas` in [output_mod.f90:439-1044](file:///home/hmonteiro/software/mocassin_HM_2026/source/output_mod.f90#L439-L1044)):**
   - With `lgCylinder2D = .true.`, `weight2D` is overridden to `1.0` (`line 581`), so `dV = V_oct`. Integrated line fluxes, $\langle T_e(X^{i+})\rangle$, and $\langle X^{i+}/X\rangle$ match 3D cylindrical averages identically.
   - With `lgCylinder2D = .false.` (planar 2D), interior cells satisfy `(Pi/2) * x_i * dx_i = (Pi/4) * (R_out^2 - R_in^2)` on uniform grids (and $i=1$ uses `xCell = dx/4`), **except at outer boundary $i = n_x$** due to Bug V4 and `xCell = xAxis(nx)`.
3. **Dust Temperatures (`getDustT` in [update_mod.f90:2446-2617](file:///home/hmonteiro/software/mocassin_HM_2026/source/update_mod.f90#L2446-L2617)):**
   - Uses local `radField = JSte(cellP,:)/Pi` (and `JDif` when `lgDebug`) and computes grain-size / species averages `Tdust(0,0,cellP)` weighted purely by `grainAbun` and `grainWeight`; independent of cell volume once `JSte` is properly normalized.

---

### 2.2 Detailed Averaging & Normalization Findings

#### `[ ]` **A1 (Medium) — `dustHeatingBudget` Omits Cell Volume `dV`**
- **Location:** [`source/update_mod.f90:2497-2505`](file:///home/hmonteiro/software/mocassin_HM_2026/source/update_mod.f90#L2497-L2505) (`getDustT`), reported in [`source/iteration_mod.f90:1032-1047`](file:///home/hmonteiro/software/mocassin_HM_2026/source/iteration_mod.f90#L1032-L1047).
- **Root Cause:** `dustHeatingBudget(abFileIndex, 0)` sums `dustAbsIntegral * grainWeight(ai) * grainAbun(nspU, nS) * grid%Ndust(cellP)` per cell without multiplying by `dV = getVolume(grid, xP, yP, zP)` (and `resLineHeating` in [`update_mod.f90:2619-2660`](file:///home/hmonteiro/software/mocassin_HM_2026/source/update_mod.f90#L2619-L2660) — check `budget` update in lines `2661-2707` as well).
- **Impact:** In `2D cylinder` grids where cell volumes scale as $\propto R_i$, omitting `dV` over-weights dense inner/axial rings ($R \to 0$) in the `Dust Budgets:` diagnostic printed to `output/summary.out` (no effect on `Tdust` or photon transfer itself).
- **Recommended Fix:** Multiply cell contributions to `dustHeatingBudget` in `update_mod.f90` by `getVolume(grid, xP, yP, zP)` (and `weight2D` in planar 2D).

#### `[x]` **A2 (Critical) — Out-of-Bounds Write on `JDif` When `lg2D .and. lgDebug` (+ Dead `TwoDscaleJ`)**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): removed dead `twoDscaleJ` (`common_mod.f90`, `grid_mod.f90`, `iteration_mod.f90`) and the out-of-bounds `JDif(:, 1:nLines)` loop in `iterateMC`; added lazy allocation of `linePDF`/`linePackets` in `iterateMC` (`nLines` is set in `makeElements` after `fillGrid`); aligned `nLines` & `linePDF` indexing with `outputGas` (`hydro_mod.f90`, `emission_mod.f90`); and added discrete CDF bisection in `getNu2` (`photon_mod.f90`). Verified under `CO=debug` (`-fbounds-check`) on 14 MPI ranks (`PN75_2D_cylinder + debug` vs non-debug: $H\beta$ $-0.01\%$, `[OIII] 5008` $+0.03\%$, Analytical/MC $=1.0125$) and full `-O2` 4-model gate ($\le 0.38\%$).
- **Location:** [`source/iteration_mod.f90:93-110,692-706`](file:///home/hmonteiro/software/mocassin_HM_2026/source/iteration_mod.f90#L93-L110), plus `TwoDscaleJ` in [`source/grid_mod.f90:1584-1585,3924-3926`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L1584-L1585).
- **Root Cause:**
  ```fortran
  if (lgDebug) then
     do freq = 1, nLines
        grid(iG)%linePackets(grid(iG)%active(i,j,k),freq) &
             & = linePacketsTemp(grid(iG)%active(i,j,k),freq)
        if (lg2D) grid(iG)%JDif(grid(iG)%active(i,j,k),freq) = &
             & grid(iG)%JDif(grid(iG)%active(i,j,k),freq)/&
             & TwoDscaleJ(grid(iG)%active(i,j,k))
     end do
  end if
  ```
  `grid(iG)%JDif` is allocated as `(0:nCells, 1:nbins)` ([grid_mod.f90:1724](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L1724)), whereas the loop runs `freq = 1, nLines` (`nLines > nbins`)! In addition, `TwoDscaleJ` is allocated only once with size `grid(1)%nCells` (out-of-bounds if `iG > 1` has more cells) and never set to anything other than `1.0`.
- **Impact:** Enabling the `debug` keyword in any 2D model causes immediate out-of-bounds heap corruption / SIGSEGV in `iterateMC`.
- **Recommended Fix:** Move any 2D `JDif` scaling to the `do freq = 1, nbins` loop above (`lines 675-686`), or remove the dead `TwoDscaleJ` division completely since `getVolume` / `cell_volumes` already handle 2D volume normalization directly.

#### `[ ]` **A3 (Medium) — Spectrograph `slit` Aperture Filtering Ignores Toroidal Geometry**
- **Location:** [`source/output_mod.f90:318-333`](file:///home/hmonteiro/software/mocassin_HM_2026/source/output_mod.f90#L318-L333) (`isSlitPass`) and [`source/mocassinPlot.f90:177-190,439-458`](file:///home/hmonteiro/software/mocassin_HM_2026/source/mocassinPlot.f90#L177-L190).
- **Root Cause:** Checks Cartesian cell center `(abs(xAxis(i)) <= d1Slit/2.) .and. (abs(yAxis(j)) <= d2Slit/2.)` (for `slitAxis == 'z'`, or permutations for `'x'`, `'y'`). In `2D XZ cylinder`, `yAxis(1) = 0`, so every ring with $R_i \le d_1/2$ satisfies `|0| <= d2Slit/2` with 100% weight even though its circular circumference $2\pi R_i$ extends far outside `|y| <= d2Slit/2`.
- **Impact:** Slit-integrated line ratios and `Hbeta` fluxes in `lineFlux-slit.out` (`2D cylinder`) fail to attenuate outer rings outside the spectrograph slit width.
- **Recommended Fix:** For `lgCylinder2D`, weight each active ring `(i, 1, k)` by the exact angular fraction $\Delta\phi(R_i, Z_k; d_1, d_2, \mathrm{slitAxis}) / (2\pi)$ of the circle $(R_i \cos\phi, R_i \sin\phi, \pm Z_k)$ falling inside the slit box (e.g., if `slitAxis == 'z'`: chord intersection of $x^2 + y^2 = R_i^2$ with $|x| \le d_1/2, |y| \le d_2/2$; if `slitAxis == 'y'`, i.e., line of sight perpendicular to symmetry axis $Z$: if $|Z_k| \le d_2/2$, fraction is $\frac{2}{\pi}\arcsin(\min(1, d_1/(2 R_i)))$).

#### `[ ]` **A4 (Low) — `writeSED` / `SEDratio` Caveats**
- **Location:** [`source/output_mod.f90:2986-3092`](file:///home/hmonteiro/software/mocassin_HM_2026/source/output_mod.f90#L2986-L3092) (`writeSED`) and [`source/photon_mod.f90:1883-1967,2647-2740,3008-3103,3142-3228,3310-3395`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L1883-L1967).
- **Findings:**
  1. **Polar Inclination Binning in `lgCylinder2D` (`viewPointPhi < 0`):** Because axisymmetric geometry depends only on $\theta = \arccos(|\hat{v}_z|)$, multiplying single-hemisphere `escapedPackets` by `4.` (`lines 3050-3051`) and dividing by `SEDratio = SUM(Luminosity(1:farIR_ind, imu)) / SUM(Luminosity(1:farIR_ind, 0))` recovers the correct isotropic solid-angle normalization **when far-IR flux is nonzero**.
  2. **Gas-Only `SEDratio` NaN (`line 3040, 3078`):** If `nAngleBins > 0` in a gas-only run whose frequency grid ends short of $150\,\mu\mathrm{m}$ (`nuMin > c/(150um)`), `SUM(Luminosity(1:farIR_ind, 0)) == 0`, producing `NaN` in inclination columns of `SED.out`. *(Pre-existing in 3D too.)*
  3. **Central Cavity `ijk == 0` Omission (`line 2992`):** Unabsorbed direct stellar packets originating in a star cell inside `R_in > 0` (`active == 0`) are stored at `escapedPackets(0, freq, imu)` cleanly by `photon_mod`, **but `writeSED` (`line 2992`) only loops `if (ijk .gt. 0)`**, dropping direct unabsorbed starlight from `SED.out` whenever the star lies inside an empty central cavity (`active(starIndeces) == 0`). *(Contrast commented legacy code at line `2767`: `do i = 0, grid(iG)%nCells`.)*
  4. **Planar 2D `noSpecular` / `buffer` SED Angular Bias:** Photons escaping perpendicularly at `y = +-deltaPerp2D/2` represent artificial truncation rather than a physical 3D surface, skewing `SED(freq, imu)` and `totalEscaped`. (Only `lgCylinder2D` gives physical 3D multi-inclination SEDs.)

---

## 3. Packet Transfer Audit

### 3.1 Verified Behavior (`lgCylinder2D` + `symmetricXYZ` in `XZ`)

1. **Exact Quadratic Cylindrical Wall Stepping ([photon_mod.f90:1409-1486,2286-2315](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L1409-L1486)):**
   - Full 3D vector state `rVec = (x, y, z)` and `vHat = (vx, vy, vz)` are evolved continuously while locking `yP = 1`.
   - Distance `ds_out` to $R_{\mathrm{out},i}$ solves $(x + s v_x)^2 + (y + s v_y)^2 = R_{\mathrm{out},i}^2$ in double precision (`dble`); when $x_P > 1$ and $B_{\mathrm{cyl}} < 0$ (moving inward), `disc_in >= 0` selects the smaller positive root `ds_in` to enter ring $x_P - 1$ (`lgTargetInner = .true.`); otherwise tangent rays cleanly advance out of `xP` (`xP + 1`) without missing tangency.
2. **Azimuthal Preservation & Equatorial Reflection ([photon_mod.f90:857-860,985-988,2085-2099,3245-3266](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L3245-L3266)):**
   - Re-emitted diffuse gas/dust packets and HG-scattered packets bypass the 3D octant reflections at $x = 0$ and $y = 0$, reflecting solely across the equatorial mirror $z = z_{\min}$ when `lgSymmetricXYZ = .true.`.
3. **Continuum Optical Depth Integration ([pathIntegration_mod.f90:15-363](file:///home/hmonteiro/software/mocassin_HM_2026/source/pathIntegration_mod.f90#L15-L363)):**
   - Audited and fixed in **`67202ee`**: `initPathTauStart` and `stepPathTauCell` track $R_{\mathrm{cyl}}$ in `lgCylinder2D` and enforce `lgHasInPlaneVel` clean exits in planar 2D.

---

### 3.2 Detailed Packet Transfer Findings

#### `[x]` **T1 (Critical) — `setResLineEscapeProb` Completely Missing 2D / `lgCylinder2D` Support**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): implemented double-precision cylindrical quadratic intersections + equatorial reflection (`lgCylinder2D`) and planar 2D perpendicular specular/escape/buffer bounds inside `setResLineEscapeProb`. Validated on dusty PN150 (`resLinesTransfer 25.`, `MdMg 0.005`): 2D cylinder vs 3D $H\beta = -1.44\%$, $\langle T_e(\mathrm{H}^+)\rangle = +0.08\%$, $[\mathrm{O\,\textsc{iii}}]\,5008 = -1.31\%$, Ly$\alpha$ heating fraction $0.5139$ vs $0.5433$ (`dusty_2D_spec` $0.5066$). Regression gate (`PN75_3D`, `HII40_3D`, `HII40_2D_cylinder`, `PN150_2D_cylinder`) passed within `<= 0.6%`.
- **Location:** [`source/emission_mod.f90:1940-2542`](file:///home/hmonteiro/software/mocassin_HM_2026/source/emission_mod.f90#L1940-L2542) (`setResLineEscapeProb`), called from [`source/emission_mod.f90:150-153`](file:///home/hmonteiro/software/mocassin_HM_2026/source/emission_mod.f90#L150-L153).
- **Root Cause:** When `lgDust = .true.` and `convPercent >= resLinesTransfer` (typically enabled on dusty nebular benchmarks), MOCASSIN traces test rays across `(ndirx+1)*(ndiry+1)*(ndirz+1)` angles in `setResLineEscapeProb` to compute `tau_mu`, resonance escape probability `fEscapeResPhotons`, and dust absorption of trapped Ly$\alpha$ and heavy-element resonance lines. `setResLineEscapeProb` has **zero `lg2D` or `lgCylinder2D` logic**.
- **Impact:** In `2D XZ` (`ny = 1`), rays start at `yAxis(1) = 0`, hit `dSy = 0` (`line 2051`), increment `yP = 2 > ny` (`line 2246`), and immediately terminate at `yP > ny` (`line 2420-2427`) on step 1! Consequently `tau_mu -> 0`, `fEscapeResPhotons -> 1`, and nearly all dust heating by trapped resonance lines (and `resLinePackets`) vanishes.
- **Recommended Fix:** Implement `lgCylinder2D` cylindrical wall ray-tracing (and planar 2D perpendicular wall specular/escape closures) inside `setResLineEscapeProb`, mirroring `photon_mod.f90:1409-1486` and `2286-2315`. (Ensure `resLineHeating` in `update_mod.f90` uses the same fixed escape probabilities.)

#### `[x]` **T2 (High) — Grazing/On-Wall `abs(dSx) < 1e-10` or `abs(dSz) < 1e-10` Jumps a Full Mean Cell Width in `lgCylinder2D`**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated `pathSegment` (`photon_mod.f90`) so on-wall steps (`dSz <= 0`, `C_out >= 0 .and. B_cyl > 0`, `C_in <= 0 .and. B_cyl < -1e-6 * sqrt(R2*A)`) return zero distance cleanly for immediate index advance without depositing `JSte` or artificial optical depth, and removed the fallback mean-cell-width substitution & `1e-10` floor. Verified on a strongly stretched S-curve `(R, Z)` `HII40` `densityFile` in `2D XZ Z cylinder` (`Hb=19.3384`, `+0.12%` vs 3D; `Te=8143.7K`, `+0.09%`; `[OIII] 5008`, `+0.76%`) and 4-model regression gate ($\le 0.74\%$).
- **Location:** [`source/photon_mod.f90:1462-1511,1727-1740`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L1462-L1511) vs [`photon_mod.f90:2185`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L2185).
- **Root Cause:** Compare 3D/planar (`lines 1530-1679`) with `lgCylinder2D` (`lines 1437-1486, 1702-1715`):
  - In 3D, when `abs(dSx) < 1e-10` (packet already on a cell wall), lines `1530-1547` snap `rVec` and immediately increment/decrement `xP` **before** reaching `lines 1718-1738` (which are then only reached when `vHat == 0`, not when `ds ~ 0`).
  - In `lgCylinder2D`, `lines 1437-1486` do NOT advance `xP` or `zP` when `ds_out < 1e-10` (for example after single-precision `rVec = rVec + dS * vHat` on line `2163` lands within $10^{-11}$ of `R_out_c`, or after a dust scattering right on a boundary!) or when `abs(dSz) < 1e-10` for `vHat%z > 0` / `zP > 1`.
  - Lines `1703-1708` then overwrite zero distance with `real((xAxis(nx) - xAxis(1)) / nx)` (**an entire mean radial cell width**) instead of advancing cleanly into the next cell, attributing that large jump's `tauCell` and `JSte` to the old cell `xP` and possibly skipping across an entire shell on a stretched/logarithmic grid while `2308` increments `xP` by only `+1`.
- **Recommended Fix:** In `lgCylinder2D` (`1437-1485`), explicitly distinguish:
  1. **Moving onto/across an interior wall within `1.d-10`:** advance `xP` (`xP + 1` or `xP - 1`) / `zP` (`zP + 1` or `zP - 1`) in a loop (or snap cleanly onto the boundary) when `dSx < 1e-10` or `dSz < 1e-10` **due to `rVec` already at the target wall**, or exit cleanly if already on outer `nx`/`nz` (`rVec >= outer_wall` moving outward);
  2. **Purely parallel velocity (`A_cyl < 1e-20` or `abs(vHat%z) <= 1e-10`):** leave `dSx = 1.e35` / `dSz = 1.e35` cleanly, removing lines `1703-1708`.

#### `[x]` **T3 (High) — Planar `bufferTransfer2D` Ignores `totalLines` and Dust Physics**
- **Status:** **FIXED & VERIFIED** on `feature/2D-mode` (2026-11-26): updated `bufferTransfer2D` (`photon_mod.f90:3715-3915`) to check dust scattering (`scaOpac / opacity` with `hg`/isotropic), sample `dustPDF` when `.not. lgGas`, sample `totalLines` optically thin line escape (`lgLine = .true.`, recording `linePackets` when `lgDebug`), update `opacity` on every bounce, and return immediately at all 6 call sites in `pathSegment`. Verified on `tau1.000` (`noGas` + `buffer`), uniform pure-dust `p0tau1` sphere (`noGas` + `buffer`), `HII40_2D_buffer` ($H\beta$ $+0.26\%$, $T_e$ $-0.20\%$), `PN150_2D_buffer` ($H\beta$ $-0.20\%$, $T_e$ $+0.02\%$), and 4-model regression gate ($\le 1.48\%$).
- **Location:** [`source/photon_mod.f90:2348-2525,3715-3915`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L3715-L3915) (`bufferTransfer2D`).
- **Root Cause:** When a packet is absorbed inside the ghost perpendicular slab (`tau_int < tau_bnd`), line `3800` unconditionally executes:
  `call getNu2(grid(gP)%recPDF(act, 1:nbins), nuP_new)`
  omitting:
  1. **Line escape fraction `grid(gP)%totalLines(act)`:** In gas nebulae, a large fraction of absorbed ionizing energy escapes as collisionally excited/recombination lines (`random <= totalLines(act)` in `newPhotonPacket`, line `920`). Forcing 100% of ghost-buffer absorptions into ionizing/continuum `recPDF` keeps absorbed energy artificially trapped in ionizing continuum!
  2. **Dust scattering & thermal dust emission:** `bufferTransfer2D` never samples `scaOpac / opacity` (HG scattering vs absorption) and never uses `dustPDF` (`nuP_new` segfaults if `lgDust .and. .not. lgGas`, where `recPDF` is unallocated!).
- **Recommended Fix:** Factor the absorption/scattering/re-emission outcome (dust `probSca`, gas `totalLines` vs `recPDF`, pure-dust `dustPDF`) so `bufferTransfer2D` enforces identical radiative microphysics inside the ghost slab as `pathSegment`.

#### `[ ]` **T4 (Low) — `"diffExt"` `posDiff%z` Typo & Discrete `(xAxis(ix), 0, zAxis(iz))` Launch Positions**
- **Location:** [`source/photon_mod.f90:161-163,224-226,878-880`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L161-L163).
- **Root Cause:**
  1. Line `163`: `posDiff%z = grid(gpLoc)%zAxis(cellLoc(2))` passes `cellLoc(2)` instead of `cellLoc(3)`! (`newPhotonPacket` at `878-880` later overwrites `positionLoc` using `difSource(1..3)`, masking the bug inside `energyPacketRun`, though the caller passes a corrupt `positionIn%z = zAxis(1)` in `2D XZ`.)
  2. For both `"diffExt"` (`878-880`) and `resLinePackets` (`224-226`), packets spawn at discrete grid node `(xAxis(ix), 0, zAxis(iz))`. In `lgCylinder2D`, spawning all ring packets at the ring midpoint $x = x_i, y = 0$ rather than sampling uniformly over the annular volume ($R = \sqrt{R_{\mathrm{in}}^2 + u_1(R_{\mathrm{out}}^2 - R_{\mathrm{in}}^2)}$, $\phi = 2\pi u_2$, $Z = Z_{\mathrm{in}} + u_3(Z_{\mathrm{out}} - Z_{\mathrm{in}})$) slightly concentrates emission at exact radius $x_i$. (Sampling within the cell volume improves smoothness; check `active(ix, 1, iz)` consistency.)
- **Recommended Fix:** Fix `cellLoc(2)` -> `cellLoc(3)` at `photon_mod.f90:163`, and consider volume-uniform `(R, phi, Z)` intra-cell jittering for `diffExt` / `resLinePackets` in `lgCylinder2D`.

#### `[ ]` **T5 / T6 (Low) — Missing Support / Startup Guards for `nGrids > 1`, Off-Axis Stars, and `fluorescence_mod`**
- **Locations:**
  - **Nested Subgrids (`nGrids > 1`):** [`source/photon_mod.f90:2166-2283,3107-3230,3271-3295`](file:///home/hmonteiro/software/mocassin_HM_2026/source/photon_mod.f90#L2166-L2283) locate mother-grid `xP(1), yP(1), zP(1)` in Cartesian `xAxis, yAxis, zAxis` rather than $R_{\mathrm{cyl}}$ in `lgCylinder2D` (`yP(1) = 1`), and subgrid exit checks (`3271-3284`) test Cartesian `rVec%x`, `rVec%y`.
  - **Off-Axis / Multiple Stars (`lgMultistars` or `starPosition%x,y /= 0`):** [`source/grid_mod.f90:4071-4090`](file:///home/hmonteiro/software/mocassin_HM_2026/source/grid_mod.f90#L4071-L4090) places stars in `XZ` at `y = 0`. In `lgCylinder2D`, only on-axis sources (`x = 0, y = 0`) are physically axisymmetric point sources (an off-axis point source at $(R_0, 0, Z_0)$ represents a luminous *ring* over $\phi \in [0, 2\pi)$ unless `initPhotonPacket` randomizes $\phi$ around $R_0$, and `photon_mod.f90:690-698` only works properly on-axis when `lgSymmetricXYZ`).
  - **X-Ray Fluorescence:** [`source/fluorescence_mod.f90:423-458,825,855-875`](file:///home/hmonteiro/software/mocassin_HM_2026/source/fluorescence_mod.f90#L423-L458) implements planar 2D `minPerp2D`/`maxPerp2D` clamping and calls `getVolume(...) / 1.e45` (Wait: note `dV` unit division `/ 1.e45` at `line 825` matches `fluorescence_mod`'s cgs convention — verify unit consistency before altering), with **no `lgCylinder2D` ray-tracing branch**.
- **Recommended Fix:** Until explicitly implemented and benchmarked, add clear startup validation in `grid_mod.f90` / `mocassin.f90` rejecting incompatible combinations with `lgCylinder2D`: `nGrids > 1`, off-axis stars (`sqrt(x^2 + y^2) > 0` without ring randomization), `lgEcho`, `lgPlaneIonization`, and `lgFluorescence` (or implement cylindrical support cleanly where needed).

---

## Audit History & Verification Log

- **2026-11-25 (`67202ee`):** Audited `source/*.f90` line-by-line across volumes (`grid_mod`, `continuum_mod`, `mocassinPlot`, `fluorescence_mod`), averaging/outputs (`iteration_mod`, `output_mod`, `update_mod`), and packet/ray transfer (`photon_mod`, `emission_mod`, `pathIntegration_mod`, `mocassinWarm`, `mocassinOutput`, `set_input_mod`). Recorded **P0** as `[x]` (fixed and verified on 14 MPI cores in `67202ee`) and logged **V1–V7, A1–A4, T1–T6** as `[ ]` for remediation.
- **2026-11-26 (T1, `187c168`):** Fixed `setResLineEscapeProb` in `source/emission_mod.f90` for `lgCylinder2D` and planar 2D modes; verified against dusty 3D/2D PN150 benchmarks and 4-model regression gate.
- **2026-11-26 (A2, `7818329`):** Removed dead `twoDscaleJ` & out-of-bounds `JDif(:, 1:nLines)` loop, added `iterateMC` lazy allocation for `linePDF`/`linePackets`, aligned `nLines`/`linePDF` (`175` HI, `34` HeI), and fixed `getNu2` discrete CDF bounds; verified under `CO=debug` (`-fbounds-check`) and `-O2` 4-model regression gate.
- **2026-11-26 (V6, `7c6e580`):** Serialized `lgCylinder2D`/`lgBuffer2D`/`lgSpecular2D` in `grid3.out`, extracted `getWeight2D` for `output_mod`/`mocassinPlot`, and fixed `mocassinWarm` pre-read loop; verified `mocassinOutput` (`0.000%`), `mocassinPlot` `grid4.out`, `mocassinWarm` warm restart, and 4-model gate.
- **2026-11-26 (V2, `ce63993`):** Added positive `[0, Rnx]` radial axis & `[-Rnz, +Rnz]` axial axis support for `lgCylinder2D` without `symmetricXYZ` (automatic & sliced 3D `densityFile`), fixed `initPathTauStart` bottom `zP` check, and verified auto/dfile non-symmetric `PN75` + 4-model gate.
- **2026-11-26 (T2, `b72a89a`):** Removed mean-cell-width jump on grazing/on-wall cylindrical steps in `photon_mod.f90` (`dSx`/`dSz >= 0` exact boundary crossing & tangency guard); verified on stretched `(R, Z)` `HII40` (`+0.12%` $H\beta$ vs 3D) and 4-model gate.
- **2026-11-26 (T3, `db8f730`):** Implemented dust scattering (`scaOpac/opacity`, `hg`/isotropic), pure-dust `dustPDF` thermal emission, and gas `totalLines` line escape in `bufferTransfer2D` (`photon_mod.f90`); verified `noGas` pure-dust `buffer` tests, `HII40`/`PN150` `buffer`, and 4-model gate.
- **2026-11-26 (V1, `5ed4207`):** Added `lgCylinder2D` toroidal volume branch to `getVolumeCon` (`continuum_mod.f90`); verified `diffuseSource` `2D XZ Z cylinder` vs `3D` and 4-model gate.
- **2026-11-26 (V4, `be486c5`):** Divided planar 2D `getVolume` and planar 2D/3D `getVolumeCon` outer boundary cell widths by `factor` (`2.0` if `lgSymmetricXYZ`) and applied outer half-cell centroid in `getWeight2D`; verified exact octant box volume (`0.000011%` error), `HII40_2D_buffer`, and 4-model gate.
- **2026-11-26 (V5, `f04a3d1`):** Multiplied `dV` by `getWeight2D` inside `fillGrid` (`gasMass`/`dustMass`), `setMotherGrid`, and `setSubGrids` in `source/grid_mod.f90`; verified `totalGasMass` and `gasMass` scaling match `2D_cylinder` to `0.00%`, and 4-model gate.
- **2026-11-26 (V3):** Enforced `(XZ, Z)` orientation and `nx >= 2, nz >= 2` for `lgCylinder2D` in `set_input_mod.f90`, closed the `else if (lg2D)` validation bypass, and guarded `getVolume`/`getVolumeCon`; verified invalid-input rejections, valid `2D XY Y buffer`, and 4-model gate.
