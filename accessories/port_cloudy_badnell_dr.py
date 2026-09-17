#!/usr/bin/env python3
"""
port_cloudy_badnell_dr.py

Extracts dielectronic recombination rate coefficient fits from local Cloudy (c25.00)
data files and source code, and ports them into the Mocassin format.

Data Sources:
1. Cloudy data/badnell_dr.dat:
   - Fit coefficients (C) and excitation energies/temperatures (E in K)
   - Origin: Prof. N. R. Badnell (University of Strathclyde, UK), version 2023-05-12
   - Isoelectronic sequences from H-like through Na-like, Mg-like, etc.
   - Ground state transitions (M = 1) for Z = 2 through Z = 30.
2. Badnell (2006, ApJ 651, L73) Table 1 & 2:
   - Dedicated DR rates for Fe 3p^q (N = 12 to 18, i.e., recombining to Fe+13 through Fe+7).
3. Gu (2003/2004) fits:
   - Parameterizations for Fe+6 (N = 20) and Fe+7 (N = 19), converted from eV to Kelvin.

Formula:
  alpha_DR(T) = T^(-3/2) * sum_{i=1}^{nfit} [ C_i * exp(-E_i / T) ]  (cm^3 s^-1)
  where T is in Kelvin, C_i in cm^3 s^-1 K^(3/2), E_i in K.

Mocassin Output Format:
  Each line contains:
  elem  ion  ncore  nfit   C(1:9)   E(1:9)
  where:
    elem  : atomic number Z (2 to 30)
    ion   : recombined ion stage index (1 = neutral, 2 = singly ionized, ..., elem-1)
    ncore : number of core electrons before recombination (N = elem - ion)
    nfit  : number of active exponential terms (<= 9)
    C(1:9): 9 fitting coefficients (unused padded with 0.0)
    E(1:9): 9 fitting energies/temperatures in K (unused padded with 0.0)
"""

import os
import sys

def port_data(cloudy_dr_file, output_file):
    if not os.path.exists(cloudy_dr_file):
        raise FileNotFoundError(f"Cloudy DR data file not found: {cloudy_dr_file}")

    with open(cloudy_dr_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    z_headers = [i for i, l in enumerate(lines) if l.startswith("Z ")]
    if len(z_headers) < 2:
        raise ValueError("Could not find both Z header sections in badnell_dr.dat")

    sec1 = lines[z_headers[0] + 1 : z_headers[1]]
    sec2 = lines[z_headers[1] + 1 :]

    data = {}

    # Parse Part 1: C coefficients
    for l in sec1:
        p = l.split()
        z, n, m = int(p[0]), int(p[1]), int(p[2])
        # Only consider Z <= 30 and ground state M = 1
        if m == 1 and z <= 30:
            c_vals = [float(x) for x in p[4:]]
            data[(z, n)] = {"c": c_vals}

    # Parse Part 2: E energies
    for l in sec2:
        p = l.split()
        z, n, m = int(p[0]), int(p[1]), int(p[2])
        if m == 1 and z <= 30:
            e_vals = [float(x) for x in p[4:]]
            if (z, n) in data:
                data[(z, n)]["e"] = e_vals

    # Fe 3p^q from Badnell (2006, ApJ 651, L73) Table 1 and 2
    cFe_q = [
        [5.636e-4, 7.390e-3, 3.635e-2, 1.693e-1, 3.315e-2, 2.288e-1, 7.316e-2],                 # N=12 (Fe 3p^2)
        [1.090e-3, 7.801e-3, 1.132e-2, 4.740e-2, 1.990e-1, 3.379e-2, 1.140e-1, 1.250e-1],       # N=13 (Fe 3p^3)
        [3.266e-3, 7.637e-3, 1.005e-2, 2.527e-2, 6.389e-2, 1.564e-1],                           # N=14 (Fe 3p^4)
        [1.074e-3, 6.080e-3, 1.887e-2, 2.540e-2, 7.580e-2, 2.773e-1],                           # N=15 (Fe 3p^5)
        [9.073e-4, 3.777e-3, 1.027e-2, 3.321e-2, 8.529e-2, 2.778e-1],                           # N=16 (Fe 3p^6)
        [5.335e-4, 1.827e-3, 4.851e-3, 2.710e-2, 8.226e-2, 3.147e-1],                           # N=17 (Fe 3p^6 4s)
        [7.421e-4, 2.526e-3, 4.605e-3, 1.489e-2, 5.891e-2, 2.318e-1]                            # N=18 (Fe 3p^6 4s^2)
    ]
    EFe_q = [
        [3.628e3, 2.432e4, 1.226e5, 4.351e5, 1.411e6, 6.589e6, 1.030e7],
        [1.246e3, 1.063e4, 4.719e4, 1.952e5, 5.637e5, 2.248e6, 7.202e6, 3.999e9],
        [1.242e3, 1.001e4, 4.466e4, 1.497e5, 3.919e5, 6.853e5],
        [1.387e3, 1.048e4, 3.955e4, 1.461e5, 4.010e5, 7.208e5],
        [1.525e3, 1.071e4, 4.033e4, 1.564e5, 4.196e5, 7.580e5],
        [2.032e3, 1.018e4, 4.638e4, 1.698e5, 4.499e5, 7.880e5],
        [3.468e3, 1.353e4, 3.690e4, 1.957e5, 4.630e5, 8.202e5]
    ]
    for idx, n in enumerate(range(12, 19)):
        data[(26, n)] = {"c": cFe_q[idx], "e": EFe_q[idx]}

    # Fe+6 (N=20) and Fe+7 (N=19) from Gu (2003/2004), converted from eV to K
    EV_TO_K = 11604.51812
    c_gu_6 = [2.50507e-11 * (EV_TO_K**1.5), 5.60226e-11 * (EV_TO_K**1.5), 1.85001e-10 * (EV_TO_K**1.5), 3.57495e-9 * (EV_TO_K**1.5), 1.66321e-7 * (EV_TO_K**1.5)]
    e_gu_6 = [8.30501e-2 * EV_TO_K, 8.52897e-1 * EV_TO_K, 3.40225e0 * EV_TO_K, 2.23053e1 * EV_TO_K, 6.80367e1 * EV_TO_K]
    data[(26, 20)] = {"c": c_gu_6, "e": e_gu_6}

    c_gu_7 = [9.19610e-11 * (EV_TO_K**1.5), 2.92460e-10 * (EV_TO_K**1.5), 1.02120e-9 * (EV_TO_K**1.5), 1.14852e-8 * (EV_TO_K**1.5), 3.25418e-7 * (EV_TO_K**1.5)]
    e_gu_7 = [1.44392e-1 * EV_TO_K, 9.23999e-1 * EV_TO_K, 5.45498e0 * EV_TO_K, 2.04301e1 * EV_TO_K, 7.06112e1 * EV_TO_K]
    data[(26, 19)] = {"c": c_gu_7, "e": e_gu_7}

    # Sort transitions by atomic number Z and recombined ion stage
    sorted_keys = sorted(data.keys(), key=lambda x: (x[0], x[0] - x[1]))

    out_lines = []
    for k in sorted_keys:
        z, ncore = k
        ion = z - ncore  # Recombined ion stage in Mocassin (ion+1 -> ion)
        c_vals = data[k]["c"]
        e_vals = data[k]["e"]
        nfit = len(c_vals)

        c_pad = c_vals + [0.0] * (9 - nfit)
        e_pad = e_vals + [0.0] * (9 - nfit)

        c_str = " ".join(f"{x:13.6E}" for x in c_pad)
        e_str = " ".join(f"{x:13.6E}" for x in e_pad)
        out_lines.append(f"{z:3d} {ion:3d} {ncore:3d} {nfit:2d}  {c_str}  {e_str}")

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, "w") as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"Successfully ported {len(out_lines)} transitions to {output_file}")


if __name__ == "__main__":
    cloudy_file = "/home/hmonteiro/software/cloudy-c25.00/data/badnell_dr.dat"
    target_file = "/home/hmonteiro/software/mocassin_HM_2026/data/badnell_dr.dat"
    port_data(cloudy_file, target_file)
