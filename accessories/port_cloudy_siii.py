#!/usr/bin/env python3
"""
port_cloudy_siii.py

Extracts [S III] atomic data from local Cloudy (c25.00) Stout database
and ports them into the MOCASSIN data/siii.dat format.

Data Sources:
1. Cloudy data/stout/s/s_3/s_3.coll:
   - Electron impact collision strengths from Hudson, C.E., Ramsbottom, C.A., Scott, M.P. 2012, ApJ, 750, 65
   - 17 temperatures from 1,000 K to 1.0e10 K
2. Cloudy data/stout/s/s_3/s_3.tp:
   - Transition probabilities (A-values) from Podobedova et al. (2009) / Froese Fischer et al. (2006)
3. Cloudy data/stout/s/s_3/s_3.nrg:
   - Energy levels and statistical weights from NIST ASD / Stout

Output:
   data/siii.dat in MOCASSIN 5-level format with irats = 0.
"""

import os
import sys

def port_siii():
    cloudy_stout_dir = os.environ.get(
        "CLOUDY_STOUT_DIR",
        "/home/hmonteiro/software/cloudy-c25.00/data/stout/s/s_3"
    )

    coll_file = os.path.join(cloudy_stout_dir, "s_3.coll")
    tp_file = os.path.join(cloudy_stout_dir, "s_3.tp")
    nrg_file = os.path.join(cloudy_stout_dir, "s_3.nrg")

    if not (os.path.exists(coll_file) and os.path.exists(tp_file) and os.path.exists(nrg_file)):
        print(f"Error: Stout files not found in {cloudy_stout_dir}", file=sys.stderr)
        sys.exit(1)

    # 1. Read temperatures and collision strengths
    with open(coll_file, "r") as f:
        f.readline()  # date/version
        temps_line = f.readline().split()
        temps = [float(x) for x in temps_line[1:]]

        cs_data = {}
        for line in f:
            if line.startswith("CS\tELECTRON"):
                parts = line.split()
                low = int(parts[2])
                up = int(parts[3])
                if low <= 5 and up <= 5:
                    vals = [float(x) for x in parts[4:]]
                    cs_data[(low, up)] = vals

    # 2. Read A-values
    A_vals = {}
    with open(tp_file, "r") as f:
        for line in f:
            if line.startswith("A\t"):
                parts = line.split()
                low = int(parts[1])
                up = int(parts[2])
                val = float(parts[3])
                if low <= 5 and up <= 5:
                    A_vals[(low, up)] = val

    # 3. Read energies and statistical weights
    levels = {}
    with open(nrg_file, "r") as f:
        f.readline()  # date/version
        for line in f:
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                idx = int(parts[0])
                if idx <= 5:
                    en = float(parts[1])
                    g = int(float(parts[2]))
                    levels[idx] = (en, g)

    # 4. Construct MOCASSIN format
    lines = []
    lines.append("4")
    lines.append(" Collision strengths: Hudson et al. (2012, ApJ, 750, 65) - Cloudy Stout")
    lines.append(" Transition probabilities: Podobedova et al. (2009) / Froese Fischer et al. (2006)")
    lines.append(" Energy levels: NIST ASD / Cloudy Stout s_3.nrg")
    lines.append(" Ported from Cloudy c25.00 for MOCASSIN")
    lines.append(f"5 {len(temps)} 0")
    lines.append("1  3p2   3P0")
    lines.append("2        3P1")
    lines.append("3        3P2")
    lines.append("4        1D2")
    lines.append("5        1S0")
    for t in temps:
        lines.append(f"{t:12.4e}")
    lines.append("0")

    for low in range(1, 5):
        for up in range(low + 1, 6):
            vals = cs_data[(low, up)]
            for k, (v, t) in enumerate(zip(vals, temps)):
                if k == 0:
                    lines.append(f"   {low:2d}  {up:2d}   {v:12.4e} {t:12.4e}")
                else:
                    lines.append(f"    0   0   {v:12.4e} {t:12.4e}")

    lines.append(" 0 0 0")

    for k in range(1, 5):
        for l in range(k + 1, 6):
            a = A_vals.get((k, l), 0.0)
            lines.append(f"   {k:2d}  {l:2d}   {a:12.4e}")

    for j in range(1, 6):
        en, g = levels[j]
        lines.append(f" {j:2d} {g:2d}  {en:12.4f}")

    output_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data",
        "siii.dat"
    )

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Successfully wrote {output_path} with {len(temps)} temperatures and {len(cs_data)} transitions.")

if __name__ == "__main__":
    port_siii()
