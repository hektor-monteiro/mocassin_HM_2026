#!/usr/bin/env python3
"""
port_cloudy_badnell_rr.py

Extracts radiative recombination rate coefficient fits from local Cloudy (c25.00)
data files and source code, and ports them into the Mocassin format.

Data Sources:
1. Cloudy data/badnell_rr.dat:
   - Fit coefficients (A, B, T0, T1, C, T2)
   - Origin: Prof. N. R. Badnell (University of Strathclyde, UK), version 2023-05-11
   - Ground state transitions (M = 1) for Z = 1 through Z = 30.
2. Badnell (2006, ApJ 651, L73) Table 3:
   - Dedicated RR rates for Fe 3p^q (N = 12 to 18, i.e., recombining to Fe+13 through Fe+7).

Formula:
  alpha_RR(T) = A / [ D * (1 + D)^(1 - B') * (1 + F)^(1 + B') ]  (cm^3 s^-1)
  where:
    D = sqrt(T / T0)
    F = sqrt(T / T1)
    B' = B + C * exp(-T2 / T)
  with T in Kelvin, T0 and T1 in Kelvin, T2 in Kelvin, A in cm^3 s^-1, B and C dimensionless.

Mocassin Output Format:
  Each line contains:
  elem  ion  ncore   A   B   T0   T1   C   T2
  where:
    elem  : atomic number Z (1 to 30)
    ion   : recombined ion stage index (1 = neutral, 2 = singly ionized, ..., elem)
    ncore : number of core electrons before recombination (N = elem - ion)
    A     : fitting parameter A (cm^3 s^-1)
    B     : fitting parameter B
    T0    : fitting parameter T0 (K)
    T1    : fitting parameter T1 (K)
    C     : fitting parameter C (0.0 if omitted)
    T2    : fitting parameter T2 (K, 0.0 if omitted)
"""

import os
import shutil
import sys

def port_data(cloudy_rr_file, output_file):
    if not os.path.exists(cloudy_rr_file):
        raise FileNotFoundError(f"Cloudy RR data file not found: {cloudy_rr_file}")

    with open(cloudy_rr_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    data = {}

    for l in lines:
        p = l.split()
        if len(p) >= 8 and p[0].isdigit():
            z, n, m = int(p[0]), int(p[1]), int(p[2])
            # Only consider Z <= 30 and ground state M = 1
            if m == 1 and z <= 30:
                a = float(p[4])
                b = float(p[5])
                t0 = float(p[6])
                t1 = float(p[7])
                c = float(p[8]) if len(p) > 8 else 0.0
                t2 = float(p[9]) if len(p) > 9 else 0.0
                data[(z, n)] = (a, b, t0, t1, c, t2)

    # Dedicated Fe 3p^q RR rate coefficients from Table 3 of Badnell (2006, ApJ, 651, L73)
    # (Fe 8+ to Fe 14+, corresponding to N = 18 down to 12 core electrons)
    fe_bad06 = {
        12: (1.179e-9, 0.7096, 4.508e2, 3.393e7, 0.0154, 3.977e6),
        13: (1.050e-9, 0.6939, 4.568e2, 3.987e7, 0.0066, 5.451e5),
        14: (9.832e-10, 0.7146, 3.597e2, 3.808e7, 0.0045, 3.952e5),
        15: (8.303e-10, 0.7156, 3.531e2, 3.554e7, 0.0132, 2.951e5),
        16: (1.052e-9, 0.7370, 1.639e2, 2.924e7, 0.0224, 4.291e5),
        17: (1.338e-9, 0.7495, 7.242e1, 2.453e7, 0.0404, 4.199e5),
        18: (1.263e-9, 0.7532, 5.209e1, 2.169e7, 0.0421, 2.917e5)
    }

    for n, pars in fe_bad06.items():
        data[(26, n)] = pars

    # Sort transitions by atomic number Z and recombined ion stage
    sorted_keys = sorted(data.keys(), key=lambda x: (x[0], x[0] - x[1]))

    out_lines = []
    for k in sorted_keys:
        z, ncore = k
        ion = z - ncore  # Recombined ion stage in Mocassin (ion=1 neutral, ..., ion=z for H-like)
        a, b, t0, t1, c, t2 = data[k]
        out_lines.append(f"{z:3d} {ion:3d} {ncore:3d}  {a:13.6E} {b:10.5f} {t0:13.6E} {t1:13.6E} {c:10.5f} {t2:13.6E}")

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, "w") as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"Successfully ported {len(out_lines)} transitions to {output_file}")

    # Also copy to data.chianty-7 and data.chianty-10 if they exist
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(output_file)))
    for extra_dir in ["data.chianty-7", "data.chianty-10"]:
        target = os.path.join(base_dir, extra_dir, "badnell_rr.dat")
        if os.path.isdir(os.path.join(base_dir, extra_dir)):
            shutil.copy2(output_file, target)
            print(f"Copied to {target}")

if __name__ == "__main__":
    cloudy_file = "/home/hmonteiro/software/cloudy-c25.00/data/badnell_rr.dat"
    target_file = "/home/hmonteiro/software/mocassin_HM_2026/data/badnell_rr.dat"
    port_data(cloudy_file, target_file)
