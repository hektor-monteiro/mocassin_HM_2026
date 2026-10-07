#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mocassin_pyneb.py
-----------------
Calculate 3D integrated emission line intensities from MOCASSIN photoionization
models using PyNeb and the CHIANTI atomic database.

Key Features:
1. Full support for Collisional Excitation Lines (CELs) using CHIANTI v10
   with the exact number of atomic energy levels (NLEVS) used by MOCASSIN.
2. Full support for Optical Recombination Lines (ORLs) for H I, He I, He II,
   C II, N II, O II, and Ne II, integrated over the parent recombining ion density.
3. Automated Abundance Discrepancy Factor (ADF) analysis in planetary nebulae
   comparing (X^{2+}/H+) from ORLs and CELs for O, C, and N.
4. Flexible line list input: custom tabular format or native MOCASSIN plot.in.
5. Direct comparison against MOCASSIN lineFlux.out analytical solutions.

Author: Antigravity (Google DeepMind) for Hektor Monteiro
Date: September 2026
"""

import os
import sys
import re
import argparse
import numpy as np

try:
    import pyneb as pn
except ImportError:
    raise ImportError(
        "PyNeb is required but not installed in this Python environment. "
        "Install it with 'pip install pyneb'."
    )

# ==============================================================================
# 1. MOCASSIN CHIANTI Level Numbers (NLEVS) Configuration
# ==============================================================================

# Standard 30 astrophysical elements
ELEMENTS = [
    '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
]

ELEM2NUM = {el.upper(): i for i, el in enumerate(ELEMENTS) if el}
NUM2ELEM = {i: el for i, el in enumerate(ELEMENTS) if el}

# Roman numerals mapping
ROMAN_TO_INT = {
    'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5,
    'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10
}

# Historical and current level counts (NLEVS) for all 116 MOCASSIN ions
DEFAULT_MOCASSIN_NLEVS = {
    "Al10": 10, "Al2": 5, "Al3": 3, "Al5": 2, "Al6": 9, "Al7": 13, "Al8": 15, "Al9": 10,
    "Ar10": 2, "Ar2": 2, "Ar3": 5, "Ar4": 5, "Ar5": 5, "Ar6": 2, "Ar7": 5, "Ar8": 3,
    "C1": 6, "C2": 10, "C3": 10, "C4": 15,
    "Ca10": 3, "Ca2": 2, "Ca4": 2, "Ca5": 5, "Ca6": 5, "Ca7": 27, "Ca8": 12, "Ca9": 5,
    "Cl10": 9, "Cl2": 5, "Cl3": 5, "Cl4": 5, "Cl5": 2, "Cl9": 2,
    "Co6": 5, "Cr3": 5, "Cr8": 2, "Cr9": 5,
    "F2": 5, "F4": 6,
    "Fe10": 2, "Fe2": 142, "Fe5": 5, "Fe6": 19, "Fe7": 9, "Fe8": 2,
    "K3": 2, "K5": 5, "K6": 5, "K7": 2, "K9": 3,
    "Mg1": 5, "Mg10": 15, "Mg2": 3, "Mg4": 2, "Mg5": 9, "Mg6": 13, "Mg7": 15, "Mg8": 10, "Mg9": 10,
    "Mn10": 5, "Mn4": 5, "Mn9": 2,
    "N1": 13, "N2": 15, "N3": 10, "N4": 10, "N5": 15,
    "Na10": 17, "Na3": 2, "Na4": 9, "Na5": 13, "Na6": 15, "Na7": 10, "Na8": 10, "Na9": 15,
    "Ne10": 9, "Ne2": 2, "Ne3": 9, "Ne4": 13, "Ne5": 15, "Ne6": 10, "Ne7": 10, "Ne8": 15, "Ne9": 17,
    "O1": 7, "O2": 13, "O3": 15, "O4": 10, "O5": 10, "O6": 15, "O8": 9,
    "P10": 15, "P3": 2, "P5": 3, "P7": 2, "P8": 9, "P9": 13,
    "S10": 13, "S2": 5, "S3": 5, "S4": 12, "S5": 5, "S6": 3, "S8": 2, "S9": 5,
    "Sc5": 2,
    "Si10": 10, "Si2": 12, "Si3": 5, "Si4": 3, "Si6": 2, "Si7": 9, "Si8": 13, "Si9": 6,
    "Ti6": 2
}


def load_mocassin_nlevs_from_dir(mocassin_dir):
    """
    Dynamically read level numbers from MOCASSIN's data/fileNames.dat and data/*.dat.
    Falls back to DEFAULT_MOCASSIN_NLEVS if not found.
    """
    nlevs = dict(DEFAULT_MOCASSIN_NLEVS)
    filenames_path = os.path.join(mocassin_dir, 'data', 'fileNames.dat')
    if not os.path.exists(filenames_path):
        return nlevs

    try:
        with open(filenames_path) as f:
            for elem in range(3, 31):
                for ion in range(1, min(elem + 2, 11)):
                    line = f.readline().strip()
                    if not line:
                        continue
                    dat_path = os.path.join(mocassin_dir, line)
                    if os.path.exists(dat_path):
                        try:
                            with open(dat_path) as df:
                                ncoms = int(df.readline().split()[0])
                                for _ in range(ncoms):
                                    df.readline()
                                hdr = df.readline().split()
                                if len(hdr) >= 2:
                                    key = f"{ELEMENTS[elem]}{ion}"
                                    nlevs[key] = int(hdr[0])
                        except Exception:
                            pass
    except Exception:
        pass
    return nlevs


def setup_pyneb_chianti(chianti_dir=None):
    """
    Configure PyNeb to use CHIANTI v10.1 database.
    Searches standard locations if chianti_dir is None.
    """
    candidates = [
        chianti_dir,
        os.environ.get('XUVTOP'),
        '/home/hmonteiro/software/chianti_v10.1',
        '/home/hmonteiro/software/cloudy-c25.00/data/chianti',
        os.path.expanduser('~/software/chianti_v10.1'),
    ]

    selected_dir = None
    for cand in candidates:
        if cand and os.path.isdir(cand):
            version_file = os.path.join(cand, 'VERSION')
            masterlist = os.path.join(cand, 'masterlist', 'masterlist.ions')
            if os.path.exists(version_file) or os.path.exists(masterlist):
                selected_dir = os.path.abspath(cand)
                break

    if selected_dir is None:
        print("[WARNING] Could not find official CHIANTI database directory. "
              "PyNeb will use its default atomic datasets.")
        return False

    pn.config.set_chianti_path(selected_dir)
    try:
        pn.atomicData.addAllChianti(replace_all=True)
        print(f"[INFO] PyNeb configured with CHIANTI (Version: {pn.config.Chianti_version}) at: {selected_dir}")
        return True
    except Exception as e:
        print(f"[WARNING] Error calling addAllChianti: {e}. PyNeb will use available datasets.")
        return False


# ==============================================================================
# 2. MOCASSIN Grid Reader
# ==============================================================================

class MocassinGrid:
    """
    Parses and holds the 3D spatial grid, physical conditions, and ionic
    abundance tables produced by a MOCASSIN photoionization model run.
    """
    def __init__(self, model_dir):
        self.model_dir = os.path.abspath(model_dir)
        self.output_dir = os.path.join(self.model_dir, 'output')
        self.input_dir = os.path.join(self.model_dir, 'input')
        
        if not os.path.isdir(self.output_dir):
            raise FileNotFoundError(f"MOCASSIN output directory not found at: {self.output_dir}")

        self._read_metadata()
        self._read_abundances()
        self._read_grid_cells()

    def _read_metadata(self):
        """Read run options and metadata from grid3.out."""
        grid3_path = os.path.join(self.output_dir, 'grid3.out')
        if not os.path.exists(grid3_path):
            raise FileNotFoundError(f"Cannot find {grid3_path}")

        with open(grid3_path) as f:
            lines = [l.strip() for l in f.readlines()]

        self.nGrids = int(lines[0].split()[0])
        self.lgSymmetricXYZ = (lines[3].split()[0].upper() == 'T')
        self.nAbComponents = int(lines[15].split()[0])
        
        # Abundance files
        self.abundance_files = []
        for i in range(self.nAbComponents):
            fname = lines[16 + i].split()[0].strip('\"\'')
            if not os.path.isabs(fname):
                fname = os.path.join(self.model_dir, fname)
            self.abundance_files.append(fname)

        # nstages is located around line 42 (after emittingGrid)
        self.nstages = 7 # default
        for idx, l in enumerate(lines):
            if 'emittingGrid' in l and idx >= 35:
                self.nstages = int(lines[idx + 1].split()[0])
                break

    def _read_abundances(self):
        """Read elemental abundances for each abundance component."""
        self.elem_abun = []
        for ab_file in self.abundance_files:
            if not os.path.exists(ab_file):
                # Check in input directory
                alt = os.path.join(self.input_dir, os.path.basename(ab_file))
                if os.path.exists(alt):
                    ab_file = alt
                else:
                    raise FileNotFoundError(f"Abundance file not found: {ab_file}")

            with open(ab_file) as f:
                ab = []
                for _ in range(30):
                    line = f.readline().split()
                    if line:
                        ab.append(float(line[0]))
                    else:
                        ab.append(0.0)
            self.elem_abun.append(ab)

        self.elem_abun = np.array(self.elem_abun) # shape (nAbComponents, 30)
        # Check active elements
        self.elem_on = np.any(self.elem_abun > 1e-12, axis=0)
        self.active_elements = [e + 1 for e in range(30) if self.elem_on[e]]

    def _read_grid_cells(self):
        """Read grid0.out, grid1.out, and grid2.out."""
        g0_path = os.path.join(self.output_dir, 'grid0.out')
        g1_path = os.path.join(self.output_dir, 'grid1.out')
        g2_path = os.path.join(self.output_dir, 'grid2.out')

        with open(g0_path) as f0, open(g1_path) as f1, open(g2_path) as f2:
            ngrids = int(f0.readline().split()[0])
            nx, ny, nz, ncells, motherP, Rout = [
                float(x) if '.' in x or 'e' in x.lower() else int(x)
                for x in f0.readline().split()
            ]
            self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
            self.Rout = float(Rout)

            self.xAxis = np.array([float(f0.readline()) for _ in range(self.nx)])
            self.yAxis = np.array([float(f0.readline()) for _ in range(self.ny)])
            self.zAxis = np.array([float(f0.readline()) for _ in range(self.nz)])

            # Cell dimensions in cm
            factor = 2.0 if self.lgSymmetricXYZ else 1.0
            dx = np.zeros(self.nx)
            dx[0] = abs(self.xAxis[1] - self.xAxis[0]) / factor
            dx[-1] = abs(self.xAxis[-1] - self.xAxis[-2]) / factor
            dx[1:-1] = abs(self.xAxis[2:] - self.xAxis[:-2]) / 2.0

            dy = np.zeros(self.ny)
            dy[0] = abs(self.yAxis[1] - self.yAxis[0]) / factor
            dy[-1] = abs(self.yAxis[-1] - self.yAxis[-2]) / factor
            dy[1:-1] = abs(self.yAxis[2:] - self.yAxis[:-2]) / 2.0

            dz = np.zeros(self.nz)
            dz[0] = abs(self.zAxis[1] - self.zAxis[0]) / factor
            dz[-1] = abs(self.zAxis[-1] - self.zAxis[-2]) / factor
            dz[1:-1] = abs(self.zAxis[2:] - self.zAxis[:-2]) / 2.0

            dV_3d = np.zeros((self.nx, self.ny, self.nz))
            for i in range(self.nx):
                for j in range(self.ny):
                    for k in range(self.nz):
                        dV_3d[i, j, k] = dx[i] * dy[j] * dz[k]

            active_cells = []
            for i in range(self.nx):
                for j in range(self.ny):
                    for k in range(self.nz):
                        l0 = f0.readline().split()
                        cellP = int(l0[0])
                        
                        l1 = f1.readline().split()
                        te = float(l1[0])
                        ne = float(l1[1])
                        hden = float(l1[2])
                        abIdx = int(l1[3]) if len(l1) > 3 else 1

                        ion_dict = {}
                        for elem in self.active_elements:
                            max_ion = min(elem + 1, self.nstages)
                            l2 = f2.readline().split()
                            ion_dict[elem] = [float(x) for x in l2[:max_ion]]

                        # Only retain active, non-zero temperature cells
                        if cellP > 0 and te > 0:
                            active_cells.append({
                                'Te': te,
                                'Ne': ne,
                                'Hden': hden,
                                'abIdx': abIdx - 1, # 0-based
                                'dV': dV_3d[i, j, k],
                                'ionDen': ion_dict
                            })

        self.n_active = len(active_cells)
        if self.n_active == 0:
            raise ValueError("No active ionized cells found in MOCASSIN grid output!")

        self.Te = np.array([c['Te'] for c in active_cells])
        self.Ne = np.array([c['Ne'] for c in active_cells])
        self.Hden = np.array([c['Hden'] for c in active_cells])
        self.abIdx = np.array([c['abIdx'] for c in active_cells])
        self.dV = np.array([c['dV'] for c in active_cells])

        # Precompute ionic densities: n(X^{stage}) = Hden * (X/H) * [X^{stage}/X] (cm^-3)
        self.ion_densities = {}
        for elem in self.active_elements:
            elem_ab = self.elem_abun[self.abIdx, elem - 1]
            max_ion = min(elem + 1, self.nstages)
            for stage in range(1, max_ion + 1):
                frac = np.array([c['ionDen'][elem][stage - 1] for c in active_cells])
                key = (elem, stage)
                self.ion_densities[key] = self.Hden * elem_ab * frac

    def get_ion_density(self, elem_num, stage):
        """
        Get 1D array of number densities n(X^{stage}) [cm^-3] across active cells.
        Returns zeros if element is not active or stage out of range.
        """
        key = (elem_num, stage)
        if key in self.ion_densities:
            return self.ion_densities[key]
        return np.zeros(self.n_active)


# ==============================================================================
# 3. Line Specification & Parser
# ==============================================================================

class LineSpec:
    """
    Holds the identification for an emission line:
    - elem: element symbol ('O', 'C', 'N', 'H', 'He', etc.)
    - spec: spectrum number (1 for neutral, 2 for singly ionized, 3 for doubly ionized)
    - wave: rest wavelength in Angstroms
    - line_type: 'CEL' (Collisional) or 'REC' (Recombination / ORL)
    - label: descriptive string
    """
    def __init__(self, elem, spec, wave, line_type='CEL', label=None):
        self.elem = elem.capitalize()
        self.spec = int(spec)
        self.wave = float(wave)
        self.line_type = line_type.upper()
        if self.line_type not in ('CEL', 'REC', 'ORL'):
            raise ValueError(f"Unknown line type: {line_type}. Must be 'CEL' or 'REC'.")
        if self.line_type == 'ORL':
            self.line_type = 'REC'
        
        self.elem_num = ELEM2NUM[self.elem.upper()]
        self.ion_str = f"{self.elem}{self.spec}"
        self.label = label or f"{self.elem} {self.spec} {self.wave:.2f}"

    def __repr__(self):
        return f"<LineSpec {self.label} ({self.line_type}) {self.wave:.2f}A>"


def parse_ion_string(ion_str):
    """
    Parse ion strings into (element, spectrum_number).
    Examples:
        'O3' -> ('O', 3)
        '[O III]' -> ('O', 3)
        'OIII' -> ('O', 3)
        'C2' -> ('C', 2)
        'H1' -> ('H', 1)
        'HeII' -> ('He', 2)
    """
    s = ion_str.strip().replace('[', '').replace(']', '').replace('\'', '').replace('\"', '')
    
    # Try Regex element + Roman numeral (e.g. 'O III' or 'OIII')
    m_roman = re.match(r'^([A-Za-z]{1,2})\s*([IVXivx]+)$', s)
    if m_roman:
        el = m_roman.group(1).capitalize()
        roman = m_roman.group(2).upper()
        if roman in ROMAN_TO_INT and el.upper() in ELEM2NUM:
            return el, ROMAN_TO_INT[roman]

    # Try Regex element + integer (e.g. 'O3' or 'O 3')
    m_int = re.match(r'^([A-Za-z]{1,2})\s*(\d+)$', s)
    if m_int:
        el = m_int.group(1).capitalize()
        sp = int(m_int.group(2))
        if el.upper() in ELEM2NUM:
            return el, sp

    raise ValueError(f"Could not parse ion string: '{ion_str}'")


def is_known_orl(elem, spec, wave):
    """Detect if transition is a prominent optical recombination line."""
    if elem == 'H' or (elem == 'He' and spec in (1, 2)):
        return True
    
    w = float(wave)
    # O II ORLs
    if elem == 'O' and spec == 2:
        if (4630 <= w <= 4680) or abs(w - 4089.29) < 1.0 or abs(w - 4596.18) < 1.0:
            return True
    # C II ORLs
    if elem == 'C' and spec == 2:
        if abs(w - 4267.0) < 5.0 or abs(w - 6578.0) < 5.0:
            return True
    # N II ORLs
    if elem == 'N' and spec == 2:
        if (5660 <= w <= 5690) or abs(w - 4630.5) < 2.0 or abs(w - 4041.3) < 2.0:
            return True
    # Ne II ORLs
    if elem == 'Ne' and spec == 2:
        if (4380 <= w <= 4450):
            return True
    return False


def parse_line_list(filepath):
    """
    Parse an input file containing emission lines.
    Supports:
    - Tabular line list: Ion  Wavelength  [Type]  [Label]
    - Native MOCASSIN plot.in format: 'line <id> <wave1> [<wave2>]'
    """
    with open(filepath) as f:
        raw_lines = f.readlines()

    line_specs = []
    
    # Check if this is a plot.in file
    is_plot_in = any(l.strip().startswith('line ') for l in raw_lines)
    
    if is_plot_in:
        line_specs = parse_plot_in_file(raw_lines)
    else:
        for idx, line in enumerate(raw_lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            
            ion_str = parts[0]
            try:
                elem, spec = parse_ion_string(ion_str)
            except Exception as e:
                print(f"[WARNING] Line {idx+1}: skipping invalid ion '{ion_str}': {e}")
                continue

            wave_str = parts[1].lower()
            if wave_str.endswith('m') or wave_str.endswith('um'):
                wave = float(wave_str.rstrip('um')) * 1e4 # microns to Angstroms
            else:
                wave = float(wave_str)

            # Check line type
            line_type = 'CEL'
            label = None
            if len(parts) >= 3:
                p2 = parts[2].upper()
                if p2 in ('CEL', 'REC', 'ORL'):
                    line_type = p2
                    if len(parts) >= 4:
                        label = " ".join(parts[3:])
                else:
                    label = " ".join(parts[2:])
            
            if len(parts) == 2 or line_type == 'CEL':
                if is_known_orl(elem, spec, wave):
                    line_type = 'REC'

            line_specs.append(LineSpec(elem, spec, wave, line_type=line_type, label=label))

    return line_specs


def parse_plot_in_file(raw_lines):
    """
    Parses MOCASSIN plot.in line list and maps line IDs / wavelengths
    to their physical ion, transition, and type.
    """
    plot_entries = []
    for line in raw_lines:
        p = line.strip().split()
        if len(p) >= 3 and p[0] == 'line':
            line_id = int(p[1])
            wave = float(p[2])
            plot_entries.append((line_id, wave))

    # MOCASSIN iLine mapping:
    # 1..175: HI recombination lines (ilow=2..8, iup=ilow+1..30)
    # 176..209: HeI recombination lines (l=1..34)
    # 210..450: HeII recombination lines
    # Heavy forbidden lines
    
    hei_waves = [
        4471.5, 2945.1, 3188.7, 3613.6, 3888.6, 3964.7, 4026.2, 4120.8, 4387.9, 4437.5,
        4471.5, 4713.2, 4921.9, 5015.7, 5047.7, 5875.7, 6678.2, 7065.3, 7281.4, 9463.6,
        10830.3, 11969.1, 12527.5, 12784.9, 12790.5, 12968.4, 15083.7, 17002.4, 18685.3,
        18697.2, 19089.4, 19543.2, 20581.3, 21120.1
    ]

    line_specs = []
    for line_id, wave in plot_entries:
        if line_id == 1:
            line_specs.append(LineSpec('H', 1, 6562.8, 'REC', 'H-alpha 6563'))
        elif line_id == 2:
            line_specs.append(LineSpec('H', 1, 4861.3, 'REC', 'H-beta 4861'))
        elif line_id == 3:
            line_specs.append(LineSpec('H', 1, 4340.5, 'REC', 'H-gamma 4341'))
        elif 176 <= line_id <= 209:
            l = line_id - 176
            w = hei_waves[l] if l < len(hei_waves) else wave
            line_specs.append(LineSpec('He', 1, w, 'REC', f'He I {w:.1f}'))
        elif line_id == 210:
            line_specs.append(LineSpec('He', 2, 1640.4, 'REC', 'He II 1640'))
        elif line_id == 212:
            line_specs.append(LineSpec('He', 2, 4685.7, 'REC', 'He II 4686'))
        elif line_id == 222:
            line_specs.append(LineSpec('He', 2, 5411.5, 'REC', 'He II 5411'))
        else:
            # Heavy ion CEL mapping based on wavelength and line ID
            elem, spec = map_heavy_cel_line(line_id, wave)
            label = f"[{elem} {spec}] {wave:.0f}"
            line_specs.append(LineSpec(elem, spec, wave, 'CEL', label))

    return line_specs


def map_heavy_cel_line(line_id, wave):
    """Identify element and spectrum number from MOCASSIN line ID and wavelength."""
    w = round(wave)
    # C lines
    if w in (8727, 8729, 8730): return ('C', 1)
    if w in (1577287, 1577288): return ('C', 2)
    if w in (1908, 1909, 1906, 1907): return ('C', 3)
    # N lines
    if w in (5198, 5199, 5200, 5201): return ('N', 1)
    if w in (6548, 6549, 6583, 6584, 6585, 5755, 5756, 1218026): return ('N', 2)
    if w in (573394, 573395): return ('N', 3)
    # O lines
    if w in (6300, 6302, 5577, 5579, 631711, 1455604): return ('O', 1)
    if w in (3726, 3727, 3729, 7320, 7321, 7322, 7330, 7331, 7332): return ('O', 2)
    if w in (5007, 5008, 4959, 4363, 4364, 518134, 883392): return ('O', 3)
    # Ne lines
    if w in (3869, 3968): return ('Ne', 3)
    # S lines
    if w in (6716, 6718, 6731, 6732): return ('S', 2)
    if w in (6312, 186811, 9069, 9072): return ('S', 3)
    # Cl lines
    if w in (5518, 5519, 5538, 5539): return ('Cl', 3)
    # Ar lines
    if w in (7136, 7137): return ('Ar', 3)
    if w in (4711, 4712, 4740, 4741): return ('Ar', 4)

    # Fallback to O3
    return ('O', 3)


# ==============================================================================
# 4. Integrated Line Intensity Calculator
# ==============================================================================

class MocassinPyNebCalculator:
    """
    Engine to integrate line emissivities across the 3D MOCASSIN model grid.
    """
    def __init__(self, model_dir, chianti_dir=None, mocassin_dir=None):
        self.grid = MocassinGrid(model_dir)
        
        # Load NLEVS matching MOCASSIN data files
        if mocassin_dir and os.path.isdir(mocassin_dir):
            self.nlevs = load_mocassin_nlevs_from_dir(mocassin_dir)
        else:
            self.nlevs = dict(DEFAULT_MOCASSIN_NLEVS)

        # Setup CHIANTI in PyNeb
        setup_pyneb_chianti(chianti_dir)

        self._atom_cache = {}
        self._rec_cache = {}

    def get_atom(self, elem, spec):
        """Get or instantiate cached pn.Atom with MOCASSIN level number."""
        key = f"{elem}{spec}"
        if key not in self._atom_cache:
            nlev = self.nlevs.get(key, 15)
            try:
                atom = pn.Atom(elem, spec, NLevels=nlev)
                self._atom_cache[key] = atom
            except Exception as e:
                print(f"[WARNING] Failed to initialize Atom({elem}, {spec}, NLevels={nlev}): {e}")
                self._atom_cache[key] = None
        return self._atom_cache[key]

    def get_rec_atom(self, elem, spec):
        """Get or instantiate cached pn.RecAtom."""
        key = f"{elem}{spec}"
        if key not in self._rec_cache:
            try:
                rec_atom = pn.RecAtom(elem, spec)
                self._rec_cache[key] = rec_atom
            except Exception as e:
                print(f"[WARNING] Failed to initialize RecAtom({elem}, {spec}): {e}")
                self._rec_cache[key] = None
        return self._rec_cache[key]

    def compute_line(self, line_spec):
        """
        Compute total integrated line luminosity (erg/s) across active cells.
        Returns dict with details: luminosity, mean_emissivity, transition, status.
        """
        Te = self.grid.Te
        Ne = self.grid.Ne
        dV = self.grid.dV

        res = {
            'elem': line_spec.elem,
            'spec': line_spec.spec,
            'wave': line_spec.wave,
            'type': line_spec.line_type,
            'label': line_spec.label,
            'luminosity': 0.0,
            'transition': None,
            'status': 'OK'
        }

        # ----------------------------------------------------------------------
        # Case A: Recombination Line (ORL)
        # ----------------------------------------------------------------------
        if line_spec.line_type == 'REC':
            rec_atom = self.get_rec_atom(line_spec.elem, line_spec.spec)
            if rec_atom is None:
                res['status'] = 'RecAtom Unavailable'
                return res

            # Match closest available wavelength in RecAtom
            matched_wave = line_spec.wave
            if hasattr(rec_atom, 'labels') and len(rec_atom.labels) > 0 and line_spec.elem not in ('H', 'HE'):
                min_dw = 1e9
                matched_k = None
                for k in rec_atom.labels:
                    s = str(k).strip()
                    if '_' in s:
                        continue
                    try:
                        clean_k = float(s.rstrip('+a-'))
                        dw = abs(clean_k - line_spec.wave)
                        if dw < min_dw:
                            min_dw = dw
                            matched_k = k
                    except Exception:
                        pass
                if matched_k is not None and min_dw <= max(5.0, 0.01 * line_spec.wave):
                    matched_wave = matched_k

            # Emissivity from PyNeb
            try:
                emis = rec_atom.getEmissivity(Te, Ne, wave=matched_wave, product=False)
                if emis is None:
                    res['status'] = 'Wavelength not found in RecAtom'
                    return res

                emis = np.nan_to_num(emis, nan=0.0)
            except Exception as e:
                res['status'] = f'Emissivity error: {e}'
                return res

            # Parent / recombining ion is stage = spec + 1
            parent_stage = line_spec.spec + 1
            n_parent = self.grid.get_ion_density(line_spec.elem_num, parent_stage)
            
            # Integrated line luminosity: L = sum(epsilon * Ne * n_parent * dV) [erg/s]
            L = np.sum(emis * Ne * n_parent * dV)
            res['luminosity'] = float(L)
            res['parent_ion'] = f"{line_spec.elem} {parent_stage}"
            return res

        # ----------------------------------------------------------------------
        # Case B: Collisional Excitation Line (CEL)
        # ----------------------------------------------------------------------
        atom = self.get_atom(line_spec.elem, line_spec.spec)
        if atom is None:
            res['status'] = 'Atom Unavailable'
            return res

        try:
            trans = atom.getTransition(wave=line_spec.wave)
            res['transition'] = trans
            emis = atom.getEmissivity(Te, Ne, wave=line_spec.wave, product=False)
            emis = np.nan_to_num(emis, nan=0.0)
        except Exception as e:
            res['status'] = f'Emissivity error: {e}'
            return res

        # Emitting ion is stage = spec
        emitting_stage = line_spec.spec
        n_ion = self.grid.get_ion_density(line_spec.elem_num, emitting_stage)

        # Integrated line luminosity: L = sum(epsilon * Ne * n_ion * dV) [erg/s]
        L = np.sum(emis * Ne * n_ion * dV)
        res['luminosity'] = float(L)
        res['emitting_ion'] = f"{line_spec.elem} {emitting_stage}"
        return res

    def calculate_lines(self, line_list):
        """
        Calculate luminosities for all lines in line_list and normalize
        relative to H-beta (4861.33 A).
        """
        # Always compute H-beta reference first
        hbeta_spec = LineSpec('H', 1, 4861.33, 'REC', 'H-beta')
        hbeta_res = self.compute_line(hbeta_spec)
        L_hb = hbeta_res['luminosity']

        results = []
        for line in line_list:
            r = self.compute_line(line)
            L = r['luminosity']
            
            if L_hb > 0:
                ratio = L / L_hb
                i_hb100 = ratio * 100.0
            else:
                ratio = np.nan
                i_hb100 = np.nan

            r['L_Hbeta'] = L_hb
            r['I_ratio'] = ratio
            r['I_Hb100'] = i_hb100
            results.append(r)

        return results

    def compute_adf_diagnostics(self, results):
        """
        Calculates Abundance Discrepancy Factors (ADF) between ORL and CEL lines.
        Compares (X^{2+}/H+) from recombination and collision lines.
        """
        by_label = {r['label']: r for r in results if r['status'] == 'OK'}
        diagnostics = []

        # Reference: H-beta
        h1 = self.get_rec_atom('H', 1)
        mean_Te = float(np.average(self.grid.Te, weights=self.grid.Ne * self.grid.get_ion_density(1, 2) * self.grid.dV))
        mean_Ne = float(np.average(self.grid.Ne, weights=self.grid.Ne * self.grid.get_ion_density(1, 2) * self.grid.dV))
        em_hb = float(h1.getEmissivity(mean_Te, mean_Ne, 4, 2))

        # --- Oxygen ADF (O2+ / H+) ---
        # CEL: [O III] 5007
        o3_5007 = None
        for r in results:
            if r['elem'] == 'O' and r['spec'] == 3 and abs(r['wave'] - 5007) < 5.0:
                o3_5007 = r
                break

        # ORL: O II 4649.13 (or multiplet sum)
        o2_4649 = None
        o2_v1_sum = 0.0
        for r in results:
            if r['elem'] == 'O' and r['spec'] == 2 and r['type'] == 'REC':
                if abs(r['wave'] - 4649.13) < 0.2:
                    o2_4649 = r
                if 4638 <= r['wave'] <= 4677:
                    o2_v1_sum += r['I_ratio']

        if o3_5007 and o2_4649 and o3_5007['I_ratio'] > 0:
            o3_atom = self.get_atom('O', 3)
            o2_rec = self.get_rec_atom('O', 2)
            em_5007 = float(o3_atom.getEmissivity(mean_Te, mean_Ne, wave=5007))
            em_4649 = float(o2_rec.getEmissivity(mean_Te, mean_Ne, wave=4649.13))
            
            abun_cel = (o3_5007['I_ratio']) * (em_hb / em_5007)
            abun_orl = (o2_4649['I_ratio']) * (em_hb / em_4649)
            adf = abun_orl / abun_cel if abun_cel > 0 else np.nan

            diagnostics.append({
                'ion': 'O(2+) / H(+)',
                'cel_line': '[O III] 5007',
                'orl_line': 'O II 4649.13',
                'abun_cel': abun_cel,
                'abun_orl': abun_orl,
                'adf': adf,
                'orl_v1_sum_ratio': o2_v1_sum
            })

        # --- Carbon ADF (C2+ / H+) ---
        c3_1908 = None
        c2_4267 = None
        for r in results:
            if r['elem'] == 'C' and r['spec'] == 3 and abs(r['wave'] - 1908) < 2.0:
                c3_1908 = r
            if r['elem'] == 'C' and r['spec'] == 2 and r['type'] == 'REC' and abs(r['wave'] - 4267) < 2.0:
                c2_4267 = r

        if c3_1908 and c2_4267 and c3_1908['I_ratio'] > 0:
            c3_atom = self.get_atom('C', 3)
            c2_rec = self.get_rec_atom('C', 2)
            em_1908 = float(c3_atom.getEmissivity(mean_Te, mean_Ne, wave=1908))
            em_4267 = float(c2_rec.getEmissivity(mean_Te, mean_Ne, wave=4267.0))

            abun_cel = (c3_1908['I_ratio']) * (em_hb / em_1908)
            abun_orl = (c2_4267['I_ratio']) * (em_hb / em_4267)
            adf = abun_orl / abun_cel if abun_cel > 0 else np.nan

            diagnostics.append({
                'ion': 'C(2+) / H(+)',
                'cel_line': 'C III] 1908',
                'orl_line': 'C II 4267.0',
                'abun_cel': abun_cel,
                'abun_orl': abun_orl,
                'adf': adf
            })

        # --- Nitrogen ADF (N2+ / H+) ---
        n3_57m = None
        n2_5679 = None
        for r in results:
            if r['elem'] == 'N' and r['spec'] == 3 and abs(r['wave'] - 573394) < 1000.0:
                n3_57m = r
            if r['elem'] == 'N' and r['spec'] == 2 and r['type'] == 'REC' and abs(r['wave'] - 5679.56) < 2.0:
                n2_5679 = r

        if n3_57m and n2_5679 and n3_57m['I_ratio'] > 0:
            n3_atom = self.get_atom('N', 3)
            n2_rec = self.get_rec_atom('N', 2)
            em_57m = float(n3_atom.getEmissivity(mean_Te, mean_Ne, wave=573394))
            em_5679 = float(n2_rec.getEmissivity(mean_Te, mean_Ne, wave=5679.56))

            abun_cel = (n3_57m['I_ratio']) * (em_hb / em_57m)
            abun_orl = (n2_5679['I_ratio']) * (em_hb / em_5679)
            adf = abun_orl / abun_cel if abun_cel > 0 else np.nan

            diagnostics.append({
                'ion': 'N(2+) / H(+)',
                'cel_line': '[N III] 57m',
                'orl_line': 'N II 5679.56',
                'abun_cel': abun_cel,
                'abun_orl': abun_orl,
                'adf': adf
            })

        return diagnostics


# ==============================================================================
# 5. Output and Reporting & MOCASSIN Comparison
# ==============================================================================

def load_mocassin_lineflux(model_dir):
    """
    Parse MOCASSIN's output/lineFlux.out for formal and analytical line ratios.
    Returns (hb_lum, moc_lines) where moc_lines is a list of dicts:
    {'elem': int, 'stage': int, 'wave': float, 'flux': float, 'type': str}
    """
    lineflux_path = os.path.join(model_dir, 'output', 'lineFlux.out')
    if not os.path.exists(lineflux_path):
        return None, []

    with open(lineflux_path) as f:
        content = f.read()

    # Integrated totals before individual components
    totals_section = content.split('Component:')[0] if 'Component:' in content else content
    lines = totals_section.split('\n')
    hb_lum = None
    moc_lines = []
    current_mode = None

    for line in lines:
        s = line.strip()
        if not s:
            continue
        if 'Hbeta [E36 erg/s]:' in s and hb_lum is None:
            try:
                hb_lum = float(s.split(':')[1].strip().split()[0]) * 1e36
            except Exception:
                pass
        if 'HI recombination lines:' in s:
            current_mode = 'HI'
            continue
        elif 'HeI recombination lines:' in s:
            current_mode = 'HeI'
            continue
        elif 'HeII recombination lines:' in s:
            current_mode = 'HeII'
            continue
        elif 'Heavy ions forbidden lines' in s:
            current_mode = 'HEAVY'
            continue
        elif 'Line Diagnostics' in s or 'Resonance Lines' in s or 'Totals integrated' in s or 'Line Ratios' in s:
            current_mode = None
            continue

        parts = s.split()
        if not current_mode:
            continue

        try:
            if current_mode == 'HI':
                # elem, stage, n_up, n_low, flux, idx, wave
                elem = int(parts[0])
                flux = float(parts[4])
                wave = float(parts[6])
                moc_lines.append({'elem': elem, 'stage': 1, 'wave': wave, 'flux': flux, 'type': 'REC'})
            elif current_mode == 'HeI':
                # elem, stage, idx, wave, flux, line_id
                elem = int(parts[0])
                wave = float(parts[3])
                flux = float(parts[4])
                moc_lines.append({'elem': elem, 'stage': 1, 'wave': wave, 'flux': flux, 'type': 'REC'})
            elif current_mode == 'HeII':
                # elem, stage, n_up, n_low, flux, idx, wave
                elem = int(parts[0])
                flux = float(parts[4])
                wave = float(parts[6])
                moc_lines.append({'elem': elem, 'stage': 2, 'wave': wave, 'flux': flux, 'type': 'REC'})
            elif current_mode == 'HEAVY':
                # elem, stage, low, up, wave, flux, line_id
                elem = int(parts[0])
                stage = int(parts[1])
                wave = float(parts[4])
                flux = float(parts[5])
                moc_lines.append({'elem': elem, 'stage': stage, 'wave': wave, 'flux': flux, 'type': 'CEL'})
        except (ValueError, IndexError):
            pass

    return hb_lum, moc_lines


def compare_with_mocassin(results, moc_lines):
    """
    Cross-matches PyNeb line results with MOCASSIN lineFlux.out entries.
    Populates 'moc_flux', 'moc_flux_hb100', and 'diff_pct' for matched lines.
    """
    for r in results:
        r['moc_flux'] = None
        r['moc_flux_hb100'] = None
        r['diff_pct'] = None

        if r['status'] != 'OK':
            continue

        elem_num = ELEM2NUM.get(r['elem'].upper())
        stage = r['spec']
        wave = r['wave']

        best_m = None
        best_dw = 1e9
        for ml in moc_lines:
            if ml['elem'] == elem_num and ml['stage'] == stage:
                dw = abs(ml['wave'] - wave)
                tol = max(5.0, 0.005 * wave)
                if dw <= tol and dw < best_dw:
                    best_dw = dw
                    best_m = ml

        if best_m is not None:
            r['moc_flux'] = best_m['flux']
            r['moc_flux_hb100'] = best_m['flux'] * 100.0
            if best_m['flux'] > 0:
                r['diff_pct'] = 100.0 * (r['I_ratio'] - best_m['flux']) / best_m['flux']


def print_results_table(results, adf_diag=None, moc_hb_lum=None):
    """Print results table to terminal."""
    if not results:
        print("No line results to display.")
        return

    L_hb = results[0].get('L_Hbeta', 0.0)
    has_moc = any(r.get('moc_flux') is not None for r in results)

    print("\n" + "=" * (135 if has_moc else 105))
    print("MOCASSIN - PyNeb Integrated Emission Line Intensities (CHIANTI v10)")
    print("=" * (135 if has_moc else 105))
    if moc_hb_lum and moc_hb_lum > 0:
        diff_hb = 100.0 * (L_hb - moc_hb_lum) / moc_hb_lum
        print(f"H-beta (4861.33 A) Luminosity: PyNeb = {L_hb:.5e} erg/s | MOCASSIN = {moc_hb_lum:.5e} erg/s (Diff: {diff_hb:+.2f}%)\n")
    else:
        print(f"H-beta (4861.33 A) Total Emitted Luminosity: {L_hb:.5e} erg/s (10^36 erg/s: {L_hb/1e36:.5f})\n")

    if has_moc:
        header = f"{'Ion':6s} {'Wavelength':14s} {'Type':5s} {'Luminosity [erg/s]':20s} {'PyNeb (Hb=1)':14s} {'PyNeb (Hb=100)':16s} {'MOC (Hb=1)':14s} {'Diff (%)':10s} {'Label':24s}"
        print(header)
        print("-" * 135)
        for r in results:
            ion_str = f"{r['elem']}{r['spec']}"
            wave_str = f"{r['wave']:10.2f} A" if r['wave'] < 10000 else f"{r['wave']/1e4:8.2f} um"
            if r['status'] == 'OK':
                moc_str = f"{r['moc_flux']:14.4e}" if r['moc_flux'] is not None else f"{'--':14s}"
                diff_str = f"{r['diff_pct']:+9.2f}%" if r['diff_pct'] is not None else f"{'--':10s}"
                print(f"{ion_str:6s} {wave_str:14s} {r['type']:5s} {r['luminosity']:20.4e} {r['I_ratio']:14.4e} {r['I_Hb100']:16.4f} {moc_str} {diff_str} {r['label']:24s}")
            else:
                print(f"{ion_str:6s} {wave_str:14s} {r['type']:5s} {'--':20s} {'--':14s} {'--':16s} {'--':14s} {'--':10s} {r['status']:24s}")
        print("-" * 135)
    else:
        header = f"{'Ion':6s} {'Wavelength':14s} {'Type':5s} {'Luminosity [erg/s]':22s} {'I / I(Hbeta)':16s} {'I (Hb=100)':12s} {'Label':24s}"
        print(header)
        print("-" * 105)
        for r in results:
            ion_str = f"{r['elem']}{r['spec']}"
            wave_str = f"{r['wave']:10.2f} A" if r['wave'] < 10000 else f"{r['wave']/1e4:8.2f} um"
            if r['status'] == 'OK':
                print(f"{ion_str:6s} {wave_str:14s} {r['type']:5s} {r['luminosity']:22.4e} {r['I_ratio']:16.4e} {r['I_Hb100']:12.4f} {r['label']:24s}")
            else:
                print(f"{ion_str:6s} {wave_str:14s} {r['type']:5s} {'--':22s} {'--':16s} {'--':12s} {r['status']:24s}")
        print("-" * 105)

    if adf_diag:
        print("\n" + "=" * 75)
        print("Abundance Discrepancy Factor (ADF) Planetary Nebula Diagnostics")
        print("=" * 75)
        print(f"{'Parent Ion':15s} {'CEL Line':15s} {'ORL Line':16s} {'(X2+/H+) CEL':14s} {'(X2+/H+) ORL':14s} {'ADF':8s}")
        print("-" * 75)
        for d in adf_diag:
            print(f"{d['ion']:15s} {d['cel_line']:15s} {d['orl_line']:16s} {d['abun_cel']:14.3e} {d['abun_orl']:14.3e} {d['adf']:8.2f}")
            if 'orl_v1_sum_ratio' in d and d['orl_v1_sum_ratio'] > 0:
                print(f"  -> Total O II Multiplet 1 (V1) Sum I(V1)/I(Hb) = {d['orl_v1_sum_ratio']:.4e} ({d['orl_v1_sum_ratio']*100:.4f} on Hb=100)")
        print("=" * 75)


def save_results_csv(results, outfile):
    """Save results to CSV with comparison columns if available."""
    has_moc = any(r.get('moc_flux') is not None for r in results)
    with open(outfile, 'w') as f:
        if has_moc:
            f.write("Ion,Wavelength_Angstrom,Type,Luminosity_erg_s,Intensity_Hb1,Intensity_Hb100,Mocassin_Hb1,Mocassin_Hb100,Diff_Percent,Status,Label\n")
            for r in results:
                m_hb1 = f"{r['moc_flux']:.6e}" if r.get('moc_flux') is not None else ""
                m_hb100 = f"{r['moc_flux_hb100']:.6f}" if r.get('moc_flux_hb100') is not None else ""
                diff = f"{r['diff_pct']:.2f}" if r.get('diff_pct') is not None else ""
                f.write(f"{r['elem']}{r['spec']},{r['wave']:.4f},{r['type']},{r['luminosity']:.6e},{r['I_ratio']:.6e},{r['I_Hb100']:.6f},{m_hb1},{m_hb100},{diff},{r['status']},{r['label']}\n")
        else:
            f.write("Ion,Wavelength_Angstrom,Type,Luminosity_erg_s,Intensity_Hb1,Intensity_Hb100,Status,Label\n")
            for r in results:
                f.write(f"{r['elem']}{r['spec']},{r['wave']:.4f},{r['type']},{r['luminosity']:.6e},{r['I_ratio']:.6e},{r['I_Hb100']:.6f},{r['status']},{r['label']}\n")
    print(f"[INFO] Results saved to CSV: {outfile}")


# ==============================================================================
# 6. Main CLI Interface
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate 3D integrated emission line intensities from MOCASSIN models with PyNeb & CHIANTI v10."
    )
    parser.add_argument(
        "--model-dir", "-m", required=True,
        help="Path to MOCASSIN model directory (containing output/ and input/)."
    )
    parser.add_argument(
        "--line-list", "-l", default=None,
        help="Path to input emission line list file (Ion Wavelength [Type] [Label])."
    )
    parser.add_argument(
        "--plot-in", "-p", default=None,
        help="Path to MOCASSIN plot.in file."
    )
    parser.add_argument(
        "--chianti-dir", "-c", default=None,
        help="Path to official CHIANTI database directory."
    )
    parser.add_argument(
        "--mocassin-dir", default=None,
        help="Path to MOCASSIN source directory to read exact NLEVS."
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Output file path to save results (.csv or .dat)."
    )
    parser.add_argument(
        "--adf", action="store_true", default=True,
        help="Display Abundance Discrepancy Factor (ADF) diagnostics for PNe (default: True)."
    )
    parser.add_argument(
        "--compare-mocassin", dest="compare_mocassin", action="store_true", default=None,
        help="Compare line intensities against MOCASSIN lineFlux.out (enabled by default if file exists)."
    )
    parser.add_argument(
        "--no-compare-mocassin", dest="compare_mocassin", action="store_false",
        help="Disable comparison against MOCASSIN lineFlux.out."
    )

    args = parser.parse_args()

    # Determine line list file
    line_file = args.line_list or args.plot_in
    if not line_file:
        # Check standard default in model directory
        def_plot_in = os.path.join(args.model_dir, 'input', 'plot.in')
        def_adf_dat = os.path.join(args.model_dir, 'input', 'lines_plot_and_adf.dat')
        if os.path.exists(def_adf_dat):
            line_file = def_adf_dat
        elif os.path.exists(def_plot_in):
            line_file = def_plot_in
        else:
            parser.error("Please specify an input line list with --line-list or --plot-in.")

    print(f"[INFO] Reading line list from: {line_file}")
    lines = parse_line_list(line_file)
    print(f"[INFO] Loaded {len(lines)} line definitions.")

    # Initialize calculator
    print(f"[INFO] Loading MOCASSIN model grid from: {args.model_dir}")
    calc = MocassinPyNebCalculator(
        model_dir=args.model_dir,
        chianti_dir=args.chianti_dir,
        mocassin_dir=args.mocassin_dir
    )

    # Compute lines
    print(f"[INFO] Calculating line intensities across {calc.grid.n_active} active 3D cells...")
    results = calc.calculate_lines(lines)

    # Compare with MOCASSIN lineFlux.out if available
    moc_hb_lum = None
    if args.compare_mocassin is not False:
        lineflux_path = os.path.join(args.model_dir, 'output', 'lineFlux.out')
        if os.path.exists(lineflux_path):
            print(f"[INFO] Loading MOCASSIN lineFlux.out for comparison: {lineflux_path}")
            moc_hb_lum, moc_lines = load_mocassin_lineflux(args.model_dir)
            if moc_lines:
                compare_with_mocassin(results, moc_lines)

    # ADF diagnostics
    adf_diag = calc.compute_adf_diagnostics(results) if args.adf else None

    # Print results
    print_results_table(results, adf_diag, moc_hb_lum)

    # Export if requested
    if args.output:
        save_results_csv(results, args.output)


if __name__ == '__main__':
    main()
