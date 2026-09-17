#!/usr/bin/env python3
"""
compare_hii40_results.py

Compares the benchmark HII40 results computed using:
1. The updated Badnell dielectronic recombination dataset (Cloudy c25.00)
2. The legacy dielectronic recombination dataset (Nussbaumer & Storey 1983-86 / Aldrovandi & Pequignot 1973)
"""

import os
import re
import sys

dir_updated = sys.argv[1] if len(sys.argv) > 1 else "benchmarks/gas/HII40/run_updated/output"
dir_legacy = sys.argv[2] if len(sys.argv) > 2 else "benchmarks/gas/HII40/run_legacy/output"


def parse_ionratio(filename):
    """
    Parses ionratio.out and returns dict: (elem, ion) -> ratio
    """
    ratios = {}
    if not os.path.exists(filename):
        return ratios
    with open(filename) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3 and p[0].isdigit() and p[1].isdigit():
                try:
                    elem = int(p[0])
                    ion = int(p[1])
                    val = float(p[2])
                    ratios[(elem, ion)] = val
                except ValueError:
                    pass
    return ratios

def parse_temperatures(filename):
    """
    Parses temperature.out and returns dict: (elem, ion) -> Te
    """
    temps = {}
    if not os.path.exists(filename):
        return temps
    with open(filename) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3 and p[0].isdigit() and p[1].isdigit():
                try:
                    elem = int(p[0])
                    ion = int(p[1])
                    val = float(p[2])
                    temps[(elem, ion)] = val
                except ValueError:
                    pass
    return temps

def parse_lineflux(filename):
    """
    Parses lineFlux.out and extracts Hbeta flux and line ratios.
    """
    lines = {}
    hbeta = 0.0
    if not os.path.exists(filename):
        return hbeta, lines
    with open(filename) as f:
        content = f.read()

    hb_match = re.search(r"Hbeta\s*\[E36\s*erg/s\]:\s*([0-9.E+-]+)", content)
    if hb_match:
        hbeta = float(hb_match.group(1))

    # Pattern for lines: [LineID] lambda value
    # or specific known lines
    patterns = [
        (r"\[NII\]\s+5755\s+([0-9.E+-]+)", "[N II] 5755"),
        (r"\[NII\]\s+6548\s+([0-9.E+-]+)", "[N II] 6548"),
        (r"\[NII\]\s+6584\s+([0-9.E+-]+)", "[N II] 6584"),
        (r"\[OII\]\s+3726\s+([0-9.E+-]+)", "[O II] 3726"),
        (r"\[OII\]\s+3729\s+([0-9.E+-]+)", "[O II] 3729"),
        (r"\[OII\]\s+7318,9\s+([0-9.E+-]+)", "[O II] 7320"),
        (r"\[OII\]\s+7330,0\s+([0-9.E+-]+)", "[O II] 7330"),
        (r"\[OIII\]\s+4363\s+([0-9.E+-]+)", "[O III] 4363"),
        (r"\[OIII\]\s+4959\s+([0-9.E+-]+)", "[O III] 4959"),
        (r"\[OIII\]\s+5008\s+([0-9.E+-]+)", "[O III] 5007"),
        (r"\[NeIII\]\s+3869\s+([0-9.E+-]+)", "[Ne III] 3869"),
        (r"\[NeIII\]\s+3967\s+([0-9.E+-]+)", "[Ne III] 3967"),
        (r"\[SII\]\s+4069\s+([0-9.E+-]+)", "[S II] 4069"),
        (r"\[SII\]\s+4076\s+([0-9.E+-]+)", "[S II] 4076"),
        (r"\[SII\]\s+6716\s+([0-9.E+-]+)", "[S II] 6716"),
        (r"\[SII\]\s+6731\s+([0-9.E+-]+)", "[S II] 6731"),
        (r"\[SIII\]\s+6312\s+([0-9.E+-]+)", "[S III] 6312"),
        (r"\[SIII\]\s+9069\s+([0-9.E+-]+)", "[S III] 9069"),
        (r"\[SIII\]\s+9532\s+([0-9.E+-]+)", "[S III] 9532"),
        (r"HeI\s*\n\s*4471\.50000\s+([0-9.E+-]+)", "He I 4471"),
        (r"5875\.66016\s+([0-9.E+-]+)", "He I 5876"),
        (r"6678\.16016\s+([0-9.E+-]+)", "He I 6678"),
        (r"HeII 4686\s+([0-9.E+-]+)", "He II 4686"),
    ]

    for pat, label in patterns:
        m = re.search(pat, content)
        if m:
            lines[label] = float(m.group(1))

    return hbeta, lines

# Parse data
ion_up = parse_ionratio(os.path.join(dir_updated, "ionratio.out"))
ion_leg = parse_ionratio(os.path.join(dir_legacy, "ionratio.out"))

temp_up = parse_temperatures(os.path.join(dir_updated, "temperature.out"))
temp_leg = parse_temperatures(os.path.join(dir_legacy, "temperature.out"))

hb_up, lines_up = parse_lineflux(os.path.join(dir_updated, "lineFlux.out"))
hb_leg, lines_leg = parse_lineflux(os.path.join(dir_legacy, "lineFlux.out"))

print("=" * 85)
print("BENCHMARK HII40: UPDATED (BADNELL) vs LEGACY (NUSSBAUMER & STOREY / ALDROVANDI)")
print("=" * 85)
print(f"Hbeta Flux [1e36 erg/s]:  Updated = {hb_up:.4f},  Legacy = {hb_leg:.4f}  (diff: {(hb_up - hb_leg)/hb_leg*100:+.2f}%)")
print()

# Table of Ionic Fractions
elem_names = {1: "H", 2: "He", 6: "C", 7: "N", 8: "O", 10: "Ne", 16: "S"}
roman = ["I", "II", "III", "IV", "V", "VI", "VII"]

print("-" * 85)
print("IONIC FRACTIONS: <X^i> / <H+>")
print(f"{'Species':<12} {'Legacy':<18} {'Updated (Badnell)':<20} {'Diff (%)':<15} {'Ratio (Up/Leg)':<15}")
print("-" * 85)

for elem in [1, 2, 6, 7, 8, 10, 16]:
    name = elem_names.get(elem, str(elem))
    max_stage = 5 if elem in [6, 7, 8, 10, 16] else (3 if elem == 2 else 2)
    for ion in range(1, max_stage + 1):
        spec = f"{name} {roman[ion-1]}"
        v_leg = ion_leg.get((elem, ion), 0.0)
        v_up = ion_up.get((elem, ion), 0.0)
        if v_leg > 1e-15 or v_up > 1e-15:
            diff_pct = (v_up - v_leg) / v_leg * 100 if v_leg > 0 else float("nan")
            ratio = v_up / v_leg if v_leg > 0 else float("nan")
            print(f"{spec:<12} {v_leg:<18.4e} {v_up:<20.4e} {diff_pct:+8.2f}%       {ratio:8.4f}")

print()
print("-" * 85)
print("MEAN ELECTRON TEMPERATURES: Te(ion) [K]")
print(f"{'Species':<12} {'Legacy [K]':<18} {'Updated [K]':<20} {'Diff (K)':<15} {'Diff (%)':<15}")
print("-" * 85)

for elem in [1, 2, 6, 7, 8, 10, 16]:
    name = elem_names.get(elem, str(elem))
    max_stage = 4 if elem in [6, 7, 8, 10, 16] else (2 if elem == 2 else 2)
    for ion in range(1, max_stage + 1):
        spec = f"{name} {roman[ion-1]}"
        t_l = temp_leg.get((elem, ion), 0.0)
        t_u = temp_up.get((elem, ion), 0.0)
        if t_l > 0 or t_u > 0:
            d_k = t_u - t_l
            d_pct = d_k / t_l * 100 if t_l > 0 else float("nan")
            print(f"{spec:<12} {t_l:<18.1f} {t_u:<20.1f} {d_k:+8.1f}        {d_pct:+6.2f}%")

print()
print("-" * 85)
print("EMISSION LINE INTENSITIES (Hbeta = 1.000)")
print(f"{'Line Identifier':<20} {'Legacy':<18} {'Updated':<18} {'Diff (%)':<15} {'Ratio (Up/Leg)':<15}")
print("-" * 85)

for label in lines_up:
    fl_leg = lines_leg.get(label, 0.0)
    fl_up = lines_up.get(label, 0.0)
    d_pct = (fl_up - fl_leg) / fl_leg * 100 if fl_leg > 0 else float("nan")
    rat = fl_up / fl_leg if fl_leg > 0 else float("nan")
    print(f"{label:<20} {fl_leg:<18.5f} {fl_up:<18.5f} {d_pct:+8.2f}%       {rat:8.4f}")

print("=" * 85)
