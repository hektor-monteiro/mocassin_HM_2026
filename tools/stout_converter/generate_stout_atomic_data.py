#!/usr/bin/env python3
"""
generate_stout_atomic_data.py

Automated converter to ingest Cloudy's Stout database files (.nrg, .tp, .coll)
and generate native MOCASSIN atomic data files (data/*.dat).

Adheres strictly to the Fortran 90 parser specification in source/hydro_mod.f90:
- Level labels formatted as (A20)
- Energy levels (cm^-1) and statistical weights g
- Einstein A-coefficients (s^-1) ordered as:
    do k = 1, NLEVS - 1
       do l = k + 1, NLEVS
          k   l   A_lk
- Thermally averaged collision strengths Upsilon(T) on temperature grid
- Enforces Upsilon >= 1.0e-30 to guard against premature exit on qx == 0.0
- Supports conversion of Stout RATE ELECTRON de-excitation rates to Upsilon(T)
- Supports level count caps (NLEVS <= 17) to respect MOCASSIN array bounds
"""

import argparse
import os
import re
import sys
import glob
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.interpolate import splrep, splev

# Physical constant for collision strength conversion (CGS)
# q_ul(T) = (8.6291e-6 / (g_u * sqrt(T))) * Upsilon(T)
COLL_CONST = 8.6291e-6

# Standard element names and symbols
ELEMENTS = [
    "", "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca",
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn"
]

ROMAN_TO_INT = {
    "i": 1, "ii": 2, "iii": 3, "iv": 4, "v": 5,
    "vi": 6, "vii": 7, "viii": 8, "ix": 9, "x": 10,
    "xi": 11, "xii": 12, "xiii": 13, "xiv": 14, "xv": 15,
    "xvi": 16, "xvii": 17, "xviii": 18, "xix": 19, "xx": 20
}

INT_TO_ROMAN = {v: k for k, v in ROMAN_TO_INT.items()}

# Historical MOCASSIN level counts (NLEVS) for ions
DEFAULT_NLEVS = {
    "ci.dat": 6, "cii.dat": 10, "ciii.dat": 10, "civ.dat": 15, "cv.dat": 15,
    "cvi.dat": 15, "ni.dat": 13, "nii.dat": 15, "niii.dat": 10, "niv.dat": 10,
    "nv.dat": 15, "nvi.dat": 15, "nvii.dat": 15, "oi.dat": 7, "oii.dat": 13,
    "oiii.dat": 15, "oiv.dat": 10, "ov.dat": 10, "ovi.dat": 15, "ovii.dat": 15,
    "oviii.dat": 9, "fi.dat": 5, "fii.dat": 5, "fiii.dat": 5, "fiv.dat": 5,
    "fv.dat": 10, "fvi.dat": 10, "fvii.dat": 10, "neii.dat": 2, "neiii.dat": 9,
    "neiv.dat": 13, "nev.dat": 15, "nevi.dat": 10, "nevii.dat": 10, "neviii.dat": 15,
    "neix.dat": 17, "nex.dat": 9, "naii.dat": 5, "naiii.dat": 2, "naiv.dat": 9,
    "nav.dat": 13, "navi.dat": 15, "navii.dat": 10, "naviii.dat": 10, "naix.dat": 15,
    "nax.dat": 17, "mgi.dat": 5, "mgii.dat": 3, "mgiii.dat": 5, "mgiv.dat": 2,
    "mgv.dat": 9, "mgvi.dat": 13, "mgvii.dat": 15, "mgviii.dat": 10, "mgix.dat": 10,
    "mgx.dat": 15, "alii.dat": 5, "aliii.dat": 3, "aliv.dat": 5, "alv.dat": 2,
    "alvi.dat": 9, "alvii.dat": 13, "alviii.dat": 15, "alix.dat": 10, "alx.dat": 10,
    "si1.dat": 5, "si2.dat": 12, "si3.dat": 5, "si4.dat": 3, "si5.dat": 5,
    "si6.dat": 2, "si7.dat": 9, "si8.dat": 13, "si9.dat": 6, "si10.dat": 10,
    "pi.dat": 5, "pii.dat": 5, "piii.dat": 2, "piv.dat": 5, "pv.dat": 3,
    "pvi.dat": 5, "pvii.dat": 2, "pviii.dat": 9, "pix.dat": 13, "px.dat": 15,
    "si.dat": 5, "sii.dat": 5, "siii.dat": 5, "siv.dat": 12, "sv.dat": 5,
    "svi.dat": 3, "svii.dat": 5, "sviii.dat": 2, "six.dat": 5, "sx.dat": 13,
    "cli.dat": 5, "clii.dat": 5, "cliii.dat": 5, "cliv.dat": 5, "clv.dat": 2,
    "clvi.dat": 5, "clvii.dat": 5, "clviii.dat": 5, "clix.dat": 2, "clx.dat": 9,
    "ari.dat": 5, "arii.dat": 2, "ariii.dat": 5, "ariv.dat": 5, "arv.dat": 5,
    "arvi.dat": 2, "arvii.dat": 5, "arviii.dat": 3, "arix.dat": 5, "arx.dat": 2,
    "ki.dat": 3, "kii.dat": 5, "kiii.dat": 2, "kiv.dat": 5, "kv.dat": 5,
    "kvi.dat": 5, "kvii.dat": 2, "kviii.dat": 5, "kix.dat": 3, "kx.dat": 5,
    "cai.dat": 5, "caii.dat": 2, "caiii.dat": 5, "caiv.dat": 2, "cav.dat": 5,
    "cavi.dat": 5, "cavii.dat": 17, "caviii.dat": 12, "caix.dat": 5, "cax.dat": 3,
    "sci.dat": 5, "scv.dat": 2, "tivi.dat": 2, "crii.dat": 5, "criii.dat": 5,
    "crvii.dat": 5, "crviii.dat": 2, "crix.dat": 5, "mniv.dat": 5, "mnviii.dat": 5,
    "mnix.dat": 2, "mnx.dat": 5, "feii.dat": 142, "feiii.dat": 15, "feiv.dat": 15,
    "fev.dat": 5, "fevi.dat": 17, "fevii.dat": 9, "feviii.dat": 2, "feix.dat": 15,
    "fex.dat": 2, "coi.dat": 5, "covi.dat": 5, "niI.dat": 5, "niII.dat": 5,
    "cui.dat": 5, "zni.dat": 5, "vvii.dat": 2, "v7.dat": 2
}

# The 14 legacy non-CHIANTI / Stout ions
LEGACY_STOUT_IONS = [
    "ar_2", "ar_6", "ca_4", "cl_9", "f_2", "f_4", "k_3", "k_7",
    "mg_1", "p_3", "sc_5", "ti_6", "v_7", "s_3"
]

# Standard 23-point logarithmic temperature grid (100 K to 2.51e6 K)
DEFAULT_TEMPS = 10.0 ** (2.0 + 0.2 * np.arange(23))


def clean_float_str(s: str) -> float:
    """Cleans theoretical or formatted NIST numbers, e.g. '[540497.7]', '2.0*', etc."""
    s = re.sub(r"[\[\]\(\)\?\*]", "", s)
    if "+x" in s:
        s = s.split("+x")[0]
    return float(s)


def resolve_stout_dir(stout_dir: Optional[str] = None) -> str:
    """Finds Cloudy Stout directory automatically."""
    candidates = [
        stout_dir,
        os.environ.get("CLOUDY_STOUT_DIR"),
        os.path.join(os.environ.get("CLOUDY_DATA_PATH", ""), "stout"),
        "/home/hmonteiro/software/cloudy-c25.00/data/stout",
        "/usr/local/share/cloudy/data/stout"
    ]
    for c in candidates:
        if c and os.path.isdir(c):
            return os.path.abspath(c)
    raise FileNotFoundError("Could not locate Cloudy Stout directory. Use --stout-dir to specify.")


def species_to_mocassin_filename(species: str) -> str:
    """Maps Stout species name (e.g. 'f_4', 'ar_2') to MOCASSIN filename (e.g. 'fiv.dat', 'arii.dat')."""
    species = species.lower().strip()
    if species.endswith(".dat"):
        return species
    
    parts = species.split("_")
    if len(parts) == 2:
        elem, ion_str = parts[0], parts[1]
        try:
            ion_num = int(ion_str)
        except ValueError:
            return f"{species}.dat"
        
        # Silicon special casing in MOCASSIN: si2.dat, si3.dat, etc.
        if elem == "si":
            return f"si{ion_num}.dat"
        
        # Vanadium special casing
        if elem == "v":
            if ion_num == 7:
                return "vvii.dat"
            if ion_num == 2:
                return "vii.dat"
        
        roman = INT_TO_ROMAN.get(ion_num, str(ion_num))
        return f"{elem}{roman}.dat"
    
    return f"{species}.dat"


def mocassin_filename_to_stout_species(fname: str) -> Optional[Tuple[str, str]]:
    """Maps MOCASSIN filename (e.g. 'fiv.dat') to (element, stout_species) e.g. ('f', 'f_4')."""
    base = os.path.basename(fname).lower()
    if base.endswith(".dat"):
        base = base[:-4]
        
    # Check special cases
    if base == "vvii" or base == "v7":
        return "v", "v_7"
    if base == "vii" or base == "v2":
        return "v", "v_2"
    
    # Silicon: si2, si3, ...
    m = re.match(r"^si([0-9]+)$", base)
    if m:
        return "si", f"si_{m.group(1)}"
    
    # General regex: elem + roman or arabic
    m = re.match(r"^([a-z]+?)(i{1,3}|iv|v|vi{0,3}|ix|x|[0-9]+)$", base)
    if m:
        elem = m.group(1)
        ion_str = m.group(2)
        if ion_str in ROMAN_TO_INT:
            ion_num = ROMAN_TO_INT[ion_str]
        else:
            ion_num = int(ion_str)
        return elem, f"{elem}_{ion_num}"
    
    return None


class StoutIonConverter:
    """Ingests Stout data for a single species and generates MOCASSIN format."""

    def __init__(self, stout_dir: str, elem: str, species: str, max_levels: int = 17):
        self.stout_dir = stout_dir
        self.elem = elem.lower()
        self.species = species.lower()
        self.max_levels = max_levels
        
        self.ion_dir = os.path.join(self.stout_dir, self.elem, self.species)
        self.nrg_path = os.path.join(self.ion_dir, f"{self.species}.nrg")
        self.tp_path = os.path.join(self.ion_dir, f"{self.species}.tp")
        self.coll_path = os.path.join(self.ion_dir, f"{self.species}.coll")
        
        if not os.path.isdir(self.ion_dir):
            raise FileNotFoundError(f"Stout species directory not found: {self.ion_dir}")
        if not (os.path.exists(self.nrg_path) and os.path.exists(self.tp_path) and os.path.exists(self.coll_path)):
            raise FileNotFoundError(f"Missing required Stout files in {self.ion_dir}")

        self.levels = []            # list of (index, energy, stat_weight, label)
        self.a_values = {}          # (low, up) -> sum(A)
        self.coll_data = {}         # (low, up) -> list of Upsilon(T)
        self.temps = []             # unified list of temperatures
        self.comments = []          # citations / header notes

    def parse(self, nlevs_override: Optional[int] = None, resample_temps: bool = False):
        """Parses all three Stout files and prepares data."""
        self._parse_nrg(nlevs_override)
        self._parse_tp()
        self._parse_coll(resample_temps)

    def _parse_nrg(self, nlevs_override: Optional[int] = None):
        """Parses .nrg energy level file."""
        with open(self.nrg_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        parsed = []
        for line in lines[1:]:  # skip magic date line
            line_str = line.strip()
            if not line_str or line_str.startswith("#") or line_str.startswith("*"):
                continue
            parts = line_str.split()
            if len(parts) >= 3 and parts[0].isdigit():
                try:
                    idx = int(parts[0])
                    energy = clean_float_str(parts[1])
                    g = clean_float_str(parts[2])
                    raw_label = " ".join(parts[3:]).strip("\"'")
                    parsed.append({
                        "orig_idx": idx,
                        "energy": energy,
                        "g": g,
                        "label": raw_label
                    })
                except Exception:
                    continue

        if not parsed:
            raise ValueError(f"No valid energy levels found in {self.nrg_path}")

        # Sort by energy
        parsed.sort(key=lambda x: x["energy"])

        # Determine NLEVS
        moc_fname = species_to_mocassin_filename(self.species)
        target_nlevs = nlevs_override or DEFAULT_NLEVS.get(moc_fname, min(len(parsed), 15))
        target_nlevs = min(target_nlevs, len(parsed), self.max_levels)

        # Ground level reference to 0.0
        e_ground = parsed[0]["energy"]
        self.levels = []
        for i in range(target_nlevs):
            item = parsed[i]
            # format 20-character label (A20)
            clean_lbl = re.sub(r'[\"\']', '', item['label']).strip()
            lbl_str = f"{i+1:2d}  {clean_lbl}"[:20].ljust(20)
            self.levels.append({
                "idx": i + 1,
                "orig_idx": item["orig_idx"],
                "energy": max(0.0, item["energy"] - e_ground),
                "g": item["g"],
                "label": lbl_str
            })

    def _parse_tp(self):
        """Parses .tp radiative transition probabilities."""
        with open(self.tp_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        orig_to_new = {lvl["orig_idx"]: lvl["idx"] for lvl in self.levels}
        nlevs = len(self.levels)
        
        self.a_values = {}
        for line in lines[1:]:
            line_str = line.strip()
            if not line_str or line_str.startswith("#") or line_str.startswith("*"):
                continue
            parts = line_str.split()
            if len(parts) >= 4 and parts[0] in ["A", "S", "G"]:
                data_type = parts[0]
                try:
                    lo_orig = int(parts[1])
                    up_orig = int(parts[2])
                    val = clean_float_str(parts[3])
                except Exception:
                    continue

                if lo_orig in orig_to_new and up_orig in orig_to_new:
                    lo = orig_to_new[lo_orig]
                    up = orig_to_new[up_orig]
                    if lo > up:
                        lo, up = up, lo
                    if up <= nlevs:
                        # In Stout, 'A' is Einstein A (s^-1). Sum multiple decay channels (M1 + E2)
                        self.a_values[(lo, up)] = self.a_values.get((lo, up), 0.0) + val

    def _parse_coll(self, resample_temps: bool = False):
        """Parses .coll collision strength or rate file."""
        with open(self.coll_path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        orig_to_new = {lvl["orig_idx"]: lvl["idx"] for lvl in self.levels}
        nlevs = len(self.levels)
        g_map = {lvl["idx"]: lvl["g"] for lvl in self.levels}

        # Extract citations / comments from coll file
        for line in lines:
            ls = line.strip()
            if ls.startswith("#Reference") or ls.startswith("# Reference") or ls.startswith("# Title") or ls.startswith("#References"):
                self.comments.append(ls.lstrip("#").strip())
            elif ls.startswith("# http") or ls.startswith("# doi"):
                self.comments.append(ls.lstrip("#").strip())
            elif ls.startswith("#") and any(yr in ls for yr in ["198", "199", "200", "201", "202"]):
                clean = ls.lstrip("#").strip()
                if len(clean) > 5 and clean not in self.comments:
                    self.comments.append(clean)

        current_temps = []
        raw_coll_transitions = []  # list of (lo, up, is_rate, [vals], [temps])

        for line in lines[1:]:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            if line_str.startswith("*"):
                # comments below stars
                continue

            parts = line_str.split()
            if parts[0] == "TEMP":
                current_temps = [clean_float_str(x) for x in parts[1:]]
            elif parts[0] in ["CS", "RATE"]:
                is_rate = (parts[0] == "RATE")
                collider = parts[1]
                if collider != "ELECTRON":
                    # MOCASSIN CEL equilibrium uses electron collisions
                    continue
                try:
                    lo_orig = int(parts[2])
                    up_orig = int(parts[3])
                    vals = [clean_float_str(x) for x in parts[4:]]
                except Exception:
                    continue

                if lo_orig in orig_to_new and up_orig in orig_to_new:
                    lo = orig_to_new[lo_orig]
                    up = orig_to_new[up_orig]
                    if lo > up:
                        lo, up = up, lo
                    if up <= nlevs and current_temps:
                        raw_coll_transitions.append((lo, up, is_rate, vals, list(current_temps)))

        if not raw_coll_transitions:
            raise ValueError(f"No valid electron collision transitions found for lowest {nlevs} levels of {self.species}")

        # Check if all transitions share identical temperature grid
        all_same_grid = True
        first_grid = raw_coll_transitions[0][4]
        for item in raw_coll_transitions[1:]:
            if item[4] != first_grid:
                all_same_grid = False
                break

        if not resample_temps and all_same_grid:
            self.temps = first_grid
        else:
            self.temps = list(DEFAULT_TEMPS)

        # Convert and resample collision strengths
        self.coll_data = {}
        for lo, up, is_rate, vals, t_grid in raw_coll_transitions:
            t_arr = np.array(t_grid, dtype=float)
            val_arr = np.array(vals, dtype=float)

            # Convert de-excitation RATE to Upsilon(T)
            if is_rate:
                g_u = g_map[up]
                # Upsilon(T) = q_ul(T) * g_u * sqrt(T) / COLL_CONST
                upsilon_arr = (val_arr * g_u * np.sqrt(t_arr)) / COLL_CONST
            else:
                upsilon_arr = val_arr

            # Resample onto self.temps if needed
            if len(self.temps) != len(t_grid) or np.any(np.abs(np.array(self.temps) - t_arr) > 1e-3):
                # Interpolate in log-log space
                valid = (t_arr > 0) & (upsilon_arr > 0)
                if np.sum(valid) >= 2:
                    log_t = np.log10(t_arr[valid])
                    log_u = np.log10(upsilon_arr[valid])
                    log_target = np.log10(np.array(self.temps))
                    # Linear interpolation clamped at edges
                    interp_u = np.interp(log_target, log_t, log_u)
                    ups_final = 10.0 ** interp_u
                else:
                    ups_final = np.full(len(self.temps), max(1.0e-30, np.mean(upsilon_arr)))
            else:
                ups_final = upsilon_arr

            # Strictly enforce positive values (guard against qx == 0.0 premature exit)
            ups_final = np.maximum(ups_final, 1.0e-30)
            self.coll_data[(lo, up)] = ups_final.tolist()

        # For any level pair (lo, up) with missing collision data, assign floor
        for lo in range(1, nlevs):
            for up in range(lo + 1, nlevs + 1):
                if (lo, up) not in self.coll_data:
                    self.coll_data[(lo, up)] = [1.0e-30] * len(self.temps)

    def generate_mocassin_file(self, output_path: str):
        """Writes the MOCASSIN formatted atomic data file."""
        nlevs = len(self.levels)
        ntemps = len(self.temps)

        lines = []

        # Section 1: Comment Header
        header_lines = [
            f"** Atomic data ported from Cloudy c25.00 Stout database: {self.species} **",
            f"Energy levels & weights: NIST ASD / Stout ({self.species}.nrg)",
            f"Radiative transition probabilities: Stout ({self.species}.tp)"
        ]
        if self.comments:
            for c in self.comments[:12]:
                header_lines.append(f"Ref: {c}"[:78])
        else:
            header_lines.append(f"Collision strengths: Stout ({self.species}.coll)")
        header_lines.append("Generated by MOCASSIN Stout Converter (tools/stout_converter)")

        lines.append(str(len(header_lines)))
        for h in header_lines:
            lines.append(h[:78])

        # Section 2: Header Line
        lines.append(f"{nlevs:5d}{ntemps:5d}")

        # Section 3: Level Labels (A20)
        for lvl in self.levels:
            lines.append(lvl["label"])

        # Section 4: Temperature Grid
        for t in self.temps:
            lines.append(f"{t:12.4e}")

        # Section 5: iRats flag (0 = collision strengths)
        lines.append("0")

        # Section 6: Collision Strengths Matrix
        for lo in range(1, nlevs):
            for up in range(lo + 1, nlevs + 1):
                vals = self.coll_data.get((lo, up), [1.0e-30] * ntemps)
                for k, (v, t) in enumerate(zip(vals, self.temps)):
                    if k == 0:
                        lines.append(f"   {lo:2d}  {up:2d}   {v:12.4e} {t:12.4e}")
                    else:
                        lines.append(f"    0   0   {v:12.4e} {t:12.4e}")

        lines.append(" 0 0 0")

        # Section 7: Radiative Transition Probabilities (A-values)
        for k in range(1, nlevs):
            for l in range(k + 1, nlevs + 1):
                a_val = self.a_values.get((k, l), 0.0)
                lines.append(f"   {k:2d}  {l:2d}   {a_val:12.4e}")

        # Section 8: Energy Levels & Statistical Weights
        for lvl in self.levels:
            lines.append(f" {lvl['idx']:2d} {int(round(lvl['g'])):2d}  {lvl['energy']:12.4f}")

        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")


def convert_species_list(
    stout_dir: str,
    species_list: List[str],
    output_dir: str,
    max_levels: int = 17,
    nlevs_overrides: Optional[Dict[str, int]] = None,
    resample_temps: bool = False
) -> Dict[str, str]:
    """Converts a list of Stout species into MOCASSIN .dat files."""
    results = {}
    nlevs_overrides = nlevs_overrides or {}

    for sp in species_list:
        sp_clean = sp.strip().lower()
        if not sp_clean:
            continue

        res = mocassin_filename_to_stout_species(sp_clean)
        if res:
            elem, stout_sp = res
        elif "_" in sp_clean:
            elem = sp_clean.split("_")[0]
            stout_sp = sp_clean
        else:
            print(f"[Warning] Could not parse species '{sp}'. Skipping.", file=sys.stderr)
            continue

        moc_fname = species_to_mocassin_filename(stout_sp)
        out_path = os.path.join(output_dir, moc_fname)

        try:
            conv = StoutIonConverter(stout_dir, elem, stout_sp, max_levels=max_levels)
            override = nlevs_overrides.get(moc_fname) or nlevs_overrides.get(stout_sp)
            conv.parse(nlevs_override=override, resample_temps=resample_temps)
            conv.generate_mocassin_file(out_path)
            results[stout_sp] = out_path
            print(f"[Success] Ported {stout_sp} ({len(conv.levels)} levels, {len(conv.temps)} temps) -> {out_path}")
            
            # If vvii.dat, also link/create v7.dat for backward compatibility
            if moc_fname == "vvii.dat":
                v7_path = os.path.join(output_dir, "v7.dat")
                conv.generate_mocassin_file(v7_path)
                print(f"          Also generated compatibility alias -> {v7_path}")
        except Exception as e:
            print(f"[Error] Failed to convert {stout_sp}: {e}", file=sys.stderr)

    return results


def main():
    parser = argparse.ArgumentParser(description="MOCASSIN Stout Database Converter")
    parser.add_argument("--stout-dir", type=str, default=None, help="Path to Cloudy Stout directory")
    parser.add_argument("--ions", "--species", type=str, default="legacy",
                        help="Comma-separated ions to convert (e.g. 'f_4,mg_1,ar_2', 'legacy', or 'all')")
    parser.add_argument("--output-dir", type=str, default="data.stout", help="Output directory (default: data.stout)")
    parser.add_argument("--max-levels", type=int, default=17, help="Global cap on level counts (default: 17)")
    parser.add_argument("--nlevs", type=str, default=None, help="Level overrides, e.g. 'f_4=5,mg_1=5'")
    parser.add_argument("--resample-grid", action="store_true", help="Resample onto standard 23-point grid")
    parser.add_argument("--copy-to-data", action="store_true", help="Copy generated files directly into data/")

    args = parser.parse_args()

    stout_dir = resolve_stout_dir(args.stout_dir)
    print(f"Using Stout directory: {stout_dir}")

    # Parse level overrides
    nlevs_map = {}
    if args.nlevs:
        for item in args.nlevs.split(","):
            if "=" in item:
                k, v = item.split("=")
                nlevs_map[k.strip().lower()] = int(v.strip())

    # Determine species list
    if args.ions.lower() == "legacy":
        species = LEGACY_STOUT_IONS
    elif args.ions.lower() == "all":
        coll_files = glob.glob(os.path.join(stout_dir, "*/*/*.coll"))
        species = [os.path.basename(cf).replace(".coll", "") for cf in coll_files]
    else:
        species = [s.strip() for s in args.ions.split(",") if s.strip()]

    print(f"Converting {len(species)} species to '{args.output_dir}'...")
    results = convert_species_list(
        stout_dir=stout_dir,
        species_list=species,
        output_dir=args.output_dir,
        max_levels=args.max_levels,
        nlevs_overrides=nlevs_map,
        resample_temps=args.resample_grid
    )

    print(f"\nCompleted: {len(results)}/{len(species)} species converted successfully.")

    if args.copy_to_data and results:
        import shutil
        data_dir = "data"
        print(f"Copying {len(results)} files to {data_dir}...")
        for sp, path in results.items():
            fname = os.path.basename(path)
            dst = os.path.join(data_dir, fname)
            shutil.copy2(path, dst)
            print(f"  Copied {fname} -> {dst}")


if __name__ == "__main__":
    main()
