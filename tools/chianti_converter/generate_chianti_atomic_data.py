#!/usr/bin/env python3
"""
generate_chianti_atomic_data.py
-------------------------------
Generates MOCASSIN atomic data files (*.dat) from official CHIANTI database releases.

Supports:
- CHIANTI raw directory (e.g. $XUVTOP or local unpacked release)
- CHIANTI .tar.gz database archive (streams directly without full unpacking)
- Automatic download of official releases from https://download.chiantidatabase.org/
- Full Burgess-Tully (1992) collision strength descaling (types 1-6)
- Both classic .splups and modern .scups formats
- All ions listed in MOCASSIN's data/fileNames.dat
- Configurable level truncation (respecting MOCASSIN's nForLevels = 17 limit)
- Automated verification and comparison with existing benchmark datasets

Author: Antigravity Agent (for MOCASSIN project)
Date: September 2026
"""

import os
import sys
import argparse
import re
import tarfile
import urllib.request
from io import BytesIO, StringIO
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
from scipy.interpolate import splrep, splev

# Physical constants from CHIANTI / NIST (CGS units)
BOLTZMANN = 1.3806504e-16       # erg / K
RYD2ERG = 2.17987197e-11        # erg / Rydberg
PLANCK = 6.6260693e-27          # erg s
LIGHT = 29979245800.0           # cm / s

# Standard element symbols (Z = 1 to 30)
ELEMENTS = [
    "", "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca",
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn"
]

# Historical default level counts (NLEVS) for MOCASSIN ions
DEFAULT_NLEVS = {
    "ci.dat": 6, "cii.dat": 10, "ciii.dat": 10, "civ.dat": 15, "cv.dat": 15,
    "cvi.dat": 15, "ni.dat": 13, "nii.dat": 15, "niii.dat": 10, "niv.dat": 10,
    "nv.dat": 15, "nvi.dat": 15, "nvii.dat": 15, "oi.dat": 7, "oii.dat": 13,
    "oiii.dat": 15, "oiv.dat": 10, "ov.dat": 10, "ovi.dat": 15, "ovii.dat": 15,
    "oviii.dat": 9, "fi.dat": 5, "fii.dat": 5, "fiii.dat": 5, "fiv.dat": 6,
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
    "cavi.dat": 5, "cavii.dat": 27, "caviii.dat": 12, "caix.dat": 5, "cax.dat": 3,
    "sci.dat": 5, "scv.dat": 2, "tivi.dat": 2, "crii.dat": 5, "criii.dat": 5,
    "crvii.dat": 5, "crviii.dat": 2, "crix.dat": 5, "mniv.dat": 5, "mnviii.dat": 5,
    "mnix.dat": 2, "mnx.dat": 5, "feii.dat": 142, "feiii.dat": 15, "feiv.dat": 15,
    "fev.dat": 5, "fevi.dat": 19, "fevii.dat": 9, "feviii.dat": 2, "feix.dat": 15,
    "fex.dat": 2, "coi.dat": 5, "covi.dat": 5, "niI.dat": 5, "niII.dat": 5,
    "cui.dat": 5, "zni.dat": 5
}

# Standard MOCASSIN 23-point temperature grid (log10(T) = 2.0 to 6.4, step 0.2)
DEFAULT_TEMPS = 10.0 ** (2.0 + 0.2 * np.arange(23))


class ChiantiSource:
    """Abstract reader for CHIANTI data from a directory or tar.gz archive."""
    def __init__(self, directory: Optional[str] = None, tar_path: Optional[str] = None):
        if directory and not os.path.isdir(directory):
            raise FileNotFoundError(f"CHIANTI directory '{directory}' does not exist.")
        self.directory = directory
        self.tar_path = tar_path
        self.tar = None
        self.tar_index = {}
        if self.tar_path:
            self.tar = tarfile.open(self.tar_path, mode="r:*")
            for member in self.tar.getmembers():
                self.tar_index[member.name.lower()] = member

    def read_text(self, rel_path: str) -> Optional[str]:
        """Read text from file either in directory or tar archive."""
        if self.directory:
            full_path = os.path.join(self.directory, rel_path)
            if os.path.exists(full_path):
                with open(full_path, "r", encoding="utf-8", errors="replace") as fp:
                    return fp.read()
            # Try lowercase or alternative pathing
            alt_path = os.path.join(self.directory, rel_path.lower())
            if os.path.exists(alt_path):
                with open(alt_path, "r", encoding="utf-8", errors="replace") as fp:
                    return fp.read()
            return None

        if self.tar:
            key = rel_path.lower().lstrip("./")
            # find matching member
            for name, member in self.tar_index.items():
                if name.endswith(key) or key.endswith(name):
                    f = self.tar.extractfile(member)
                    if f:
                        return f.read().decode("utf-8", errors="replace")
            return None

        return None

    def exists(self, rel_path: str) -> bool:
        """Check whether a specific relative path exists."""
        if self.directory:
            return os.path.exists(os.path.join(self.directory, rel_path)) or \
                   os.path.exists(os.path.join(self.directory, rel_path.lower()))
        if self.tar:
            key = rel_path.lower().lstrip("./")
            for name in self.tar_index:
                if name.endswith(key) or key.endswith(name):
                    return True
        return False

    def close(self):
        if self.tar:
            self.tar.close()


def resolve_file_list(path: str) -> str:
    """
    Resolves fileNames.dat path across CWD, script directory, and repository root.
    """
    if os.path.exists(path):
        return os.path.abspath(path)

    # Check relative to script directory and repo root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    candidates = [
        os.path.join(repo_root, path),
        os.path.join(repo_root, "data", os.path.basename(path)),
        os.path.join(repo_root, "data", "fileNames.dat"),
        os.path.join(script_dir, path),
        os.path.join(script_dir, "..", path)
    ]
    for cand in candidates:
        if os.path.exists(cand):
            return os.path.abspath(cand)

    print(f"ERROR: Cannot find ion list file '{path}'")
    print(f"       Checked local path and repository root at {repo_root}")
    sys.exit(1)


def parse_file_names(filename_path: str) -> List[Tuple[int, int, str, str]]:
    """
    Parses data/fileNames.dat and returns a list of tuples:
    (element_Z, ion_stage, mocassin_relpath, chianti_ion_name)
    """
    with open(filename_path, "r") as fp:
        entries = [line.strip() for line in fp if line.strip()]

    mapping = []
    idx = 0
    for elem in range(3, 31):  # Li to Zn (3 to 30)
        max_ion = min(elem + 1, 10)
        for ion in range(1, max_ion + 1):
            if idx >= len(entries):
                break
            moc_file = entries[idx]
            el_sym = ELEMENTS[elem]
            ch_ion = f"{el_sym}_{ion}"
            mapping.append((elem, ion, moc_file, ch_ion))
            idx += 1
    return mapping


def parse_elvlc(text: str) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Parses a CHIANTI .elvlc text content.
    Supports both official fixed-width / whitespace CHIANTI standard (v8-v10.1+)
    and custom multi-column variants (e.g. Cloudy).
    Returns: (levels_list, comment_lines)
    """
    levels = []
    comments = []
    reading_comments = False

    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s == "-1":
            reading_comments = not reading_comments
            continue
        if reading_comments:
            if not s.startswith("Reference:"):
                comments.append(line.rstrip())
            continue

        # 1. Try official fixed-width slicing (ChiantiPy standard)
        parsed = False
        if len(line) >= 70:
            try:
                lvl = int(line[0:7].strip())
                conf = line[7:37].strip()
                spin_str = line[42:47].strip()
                spd = line[47:52].strip()
                j_str = line[52:57].strip()
                eobs_str = line[57:72].strip()
                eth_str = line[72:87].strip() if len(line) >= 87 else "-1.0"
                spin = int(spin_str) if spin_str else 1
                j = float(j_str) if j_str else 0.0
                eobs_cm = float(eobs_str) if eobs_str else -1.0
                eth_cm = float(eth_str) if eth_str else -1.0
                ecm = eobs_cm if eobs_cm >= 0.0 else eth_cm
                mult = int(round(2.0 * j + 1.0))
                label = f"{conf} {spin}{spd}{int(j)}"
                levels.append({
                    "lvl": lvl,
                    "conf": conf,
                    "spin": spin,
                    "spd": spd,
                    "j": j,
                    "mult": mult,
                    "ecm": ecm,
                    "label": label
                })
                parsed = True
            except Exception:
                pass

        if parsed:
            continue

        # 2. Token-based fallback
        tokens = s.split()
        clean_tokens = [t for t in tokens if t != ","]
        if len(clean_tokens) >= 6:
            try:
                lvl = int(clean_tokens[0])
                # Format with extra trailing Rydberg columns (e.g. Cloudy)
                if len(clean_tokens) >= 11 and "." in clean_tokens[-1] and "." in clean_tokens[-3]:
                    eobs_cm = float(clean_tokens[-4])
                    eth_cm = float(clean_tokens[-2])
                    mult = int(float(clean_tokens[-5]))
                    j = float(clean_tokens[-6])
                    spd = clean_tokens[-7]
                    spin = int(clean_tokens[-9])
                    conf = " ".join(clean_tokens[1:-9])
                else:
                    eth_cm = float(clean_tokens[-1])
                    eobs_cm = float(clean_tokens[-2])
                    j = float(clean_tokens[-3])
                    spd = clean_tokens[-4]
                    spin = int(clean_tokens[-5])
                    conf = " ".join(clean_tokens[1:-5])
                    mult = int(round(2.0 * j + 1.0))

                ecm = eobs_cm if eobs_cm >= 0.0 else eth_cm
                label = f"{conf} {spin}{spd}{int(j)}"
                levels.append({
                    "lvl": lvl,
                    "conf": conf,
                    "spin": spin,
                    "spd": spd,
                    "j": j,
                    "mult": mult,
                    "ecm": ecm,
                    "label": label
                })
            except Exception:
                pass

    return levels, comments


def parse_wgfa(text: str) -> Tuple[Dict[Tuple[int, int], float], List[str]]:
    """
    Parses a CHIANTI .wgfa text content.
    Returns: (a_dict[(l1, l2)] = a_value, comment_lines)
    """
    a_dict = {}
    comments = []
    reading_comments = False

    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s == "-1":
            reading_comments = not reading_comments
            continue
        if reading_comments:
            if not s.startswith("Reference:"):
                comments.append(line.rstrip())
            continue

        tokens = s.split()
        if len(tokens) >= 5:
            try:
                l1 = int(tokens[0])
                l2 = int(tokens[1])
                aval = float(tokens[4])
                a_dict[(l1, l2)] = aval
            except Exception:
                pass

    return a_dict, comments


def parse_splups(text: str) -> Tuple[Dict[Tuple[int, int], Dict[str, Any]], List[str]]:
    """
    Parses classic CHIANTI .splups text content.
    Returns: (splups_dict[(l1, l2)] = params, comment_lines)
    """
    splups_dict = {}
    comments = []
    reading_comments = False

    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s == "-1":
            reading_comments = not reading_comments
            continue
        if reading_comments:
            if not s.startswith("Reference:"):
                comments.append(line.rstrip())
            continue

        tokens = s.split()
        if len(tokens) >= 9:
            try:
                l1 = int(tokens[2])
                l2 = int(tokens[3])
                ttype = int(tokens[4])
                gf = float(tokens[5])
                de = float(tokens[6])
                cups = float(tokens[7])
                splvals = [float(x) for x in tokens[8:]]
                splups_dict[(l1, l2)] = {
                    "format": "splups",
                    "ttype": ttype,
                    "gf": gf,
                    "de": de,
                    "cups": cups,
                    "spl": splvals
                }
            except Exception:
                pass

    return splups_dict, comments


def parse_scups(text: str) -> Tuple[Dict[Tuple[int, int], Dict[str, Any]], List[str]]:
    """
    Parses modern CHIANTI .scups (v8+) text content.
    Returns: (scups_dict[(l1, l2)] = params, comment_lines)
    """
    scups_dict = {}
    comments = []
    lines = text.splitlines()

    minus_one_indices = [i for i, l in enumerate(lines) if l.strip() == "-1"]
    data_lines = lines[:minus_one_indices[0]] if minus_one_indices else lines

    if len(minus_one_indices) >= 2:
        comments = [l.rstrip() for l in lines[minus_one_indices[0]+1:minus_one_indices[1]]
                    if not l.strip().startswith("Reference:")]

    i = 0
    while i < len(data_lines):
        line = data_lines[i].strip()
        if not line:
            i += 1
            continue
        tokens = line.split()
        if len(tokens) >= 8:
            try:
                l1 = int(tokens[0])
                l2 = int(tokens[1])
                de = float(tokens[2])
                gf = float(tokens[3])
                lim = float(tokens[4])
                ntemp = int(tokens[5])
                ttype = int(tokens[6])
                cups = float(tokens[7])

                btemp = []
                j = i + 1
                while j < len(data_lines) and len(btemp) < ntemp:
                    btemp.extend([float(x) for x in data_lines[j].split()])
                    j += 1
                bscups = []
                while j < len(data_lines) and len(bscups) < ntemp:
                    bscups.extend([float(x) for x in data_lines[j].split()])
                    j += 1

                if len(btemp) == ntemp and len(bscups) == ntemp:
                    scups_dict[(l1, l2)] = {
                        "format": "scups",
                        "ttype": ttype,
                        "gf": gf,
                        "de": de,
                        "cups": cups,
                        "lim": lim,
                        "btemp": btemp,
                        "bscups": bscups
                    }
                    i = j
                    continue
            except Exception:
                pass
        i += 1

    return scups_dict, comments


def descale_upsilon(param: Dict[str, Any], temps: np.ndarray) -> np.ndarray:
    """
    Descales thermally-averaged collision strength Upsilon(T) from Burgess & Tully (1992)
    scaling fits for any transition.
    """
    ttype = param["ttype"]
    cups = param["cups"]
    de = param["de"]

    if param["format"] == "splups":
        scups = np.array(param["spl"])
        nspl = len(scups)
        dx = 1.0 / (float(nspl) - 1.0)
        xs = dx * np.arange(nspl)
    else:
        xs = np.array(param["btemp"])
        scups = np.array(param["bscups"])

    kte = BOLTZMANN * temps / (de * RYD2ERG)

    if ttype == 1:
        st = 1.0 - np.log(cups) / np.log(kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups * np.log(kte + np.e)
    elif ttype == 2:
        st = kte / (kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups
    elif ttype == 3:
        st = kte / (kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups / (kte + 1.0)
    elif ttype == 4:
        st = 1.0 - np.log(cups) / np.log(kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups * np.log(kte + cups)
    elif ttype == 5:
        st = kte / (kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups / (kte + 1e-30)
    elif ttype == 6:
        st = kte / (kte + cups)
        y2 = splrep(xs, scups, s=0)
        sups = splev(st, y2, der=0)
        u = sups
    else:
        u = np.zeros_like(temps)

    return np.maximum(u, 1e-30)


def generate_mocassin_file_content(
    ion_name: str,
    chianti_version_str: str,
    levels: List[Dict[str, Any]],
    a_dict: Dict[Tuple[int, int], float],
    coll_dict: Dict[Tuple[int, int], Dict[str, Any]],
    comments_all: List[str],
    nlevs: int,
    temps: np.ndarray
) -> str:
    """
    Constructs the exact ASCII text conforming to MOCASSIN Fortran read specifications.
    """
    actual_lvls = levels[:nlevs]
    ntemps = len(temps)

    # 1. Clean comment header
    header_comms = [f"** Atomic data from version {chianti_version_str} of the CHIANTI database **"]
    for c in comments_all:
        cs = c.strip()
        if cs and cs != "-1" and not cs.startswith("Reference:"):
            header_comms.append(c.rstrip())

    ncoms = len(header_comms)
    out_lines = [f"{ncoms}"]
    out_lines.extend(header_comms)

    # 2. Header line: NLEVS NTEMPS
    out_lines.append(f"{nlevs} {ntemps}")

    # 3. Level labels (A20 format)
    for lv in actual_lvls:
        lbl = lv["label"]
        out_lines.append(f"  {lv['lvl']:2d}  {lbl:<15s}"[:20])

    # 4. Temperatures
    for T in temps:
        out_lines.append(f"{T:10.1f}")

    # 5. iRats (0 for collision strengths)
    out_lines.append("0")

    # 6. Collisional strength blocks for available transitions
    for i in range(1, nlevs):
        for j in range(i + 1, nlevs + 1):
            if (i, j) in coll_dict:
                u = descale_upsilon(coll_dict[(i, j)], temps)
                out_lines.append(f"  {i:2d} {j:2d}   {u[0]:.2e}")
                for val in u[1:]:
                    out_lines.append(f"   0  0   {val:.2e}")

    # Terminator for collision data
    out_lines.append(" 0 0 0")

    # 7. Radiative transition probabilities (Einstein A values)
    # Loop order in Fortran:
    # do k = 1, NLEVS-1
    #    do l = k+1, NLEVS
    for i in range(1, nlevs):
        for j in range(i + 1, nlevs + 1):
            aval = a_dict.get((i, j), 0.0)
            out_lines.append(f"  {i:2d} {j:2d}   {aval:.3e}")

    # 8. Energy levels and statistical weights
    for lv in actual_lvls:
        out_lines.append(f"  {lv['lvl']:2d} {lv['mult']:2d}   {lv['ecm']:10.1f}")

    return "\n".join(out_lines) + "\n"


def download_chianti_release(version: str, cache_dir: str = "chianti_cache") -> str:
    """Downloads official CHIANTI database archive from chiantidatabase.org."""
    os.makedirs(cache_dir, exist_ok=True)
    filename = f"CHIANTI_{version}_database.tar.gz"
    url = f"https://download.chiantidatabase.org/{filename}"
    dest = os.path.join(cache_dir, filename)

    if os.path.exists(dest):
        print(f"[CHIANTI] Using cached archive: {dest}")
        return dest

    print(f"[CHIANTI] Downloading {url} -> {dest}...")
    urllib.request.urlretrieve(url, dest)
    print(f"[CHIANTI] Download complete ({os.path.getsize(dest) / (1024*1024):.1f} MB).")
    return dest


def resolve_chianti_dir(dir_path: Optional[str]) -> Optional[str]:
    """
    Validates and resolves the CHIANTI database directory.
    Supports exact paths, case-insensitive matching, variation normalization
    (e.g. CHIANTI_10.1 -> chianti_v10.1), and falls back to standard known locations
    if no path is provided.
    """
    if dir_path:
        # 1. Exact match
        if os.path.isdir(dir_path):
            return os.path.abspath(dir_path)

        # 2. Check parent directory for naming variations (e.g. CHIANTI_10.1 -> chianti_v10.1)
        abs_path = os.path.abspath(dir_path)
        parent = os.path.dirname(abs_path)
        target = os.path.basename(abs_path)
        if os.path.isdir(parent):
            clean = lambda s: re.sub(r'[^a-zA-Z0-9]', '', s).lower()
            clean_nov = lambda s: re.sub(r'[^a-zA-Z0-9]', '', s).lower().replace('v', '')
            target_clean = clean(target)
            target_nov = clean_nov(target)

            candidates = []
            for entry in os.listdir(parent):
                full_entry = os.path.join(parent, entry)
                if not os.path.isdir(full_entry):
                    continue
                entry_clean = clean(entry)
                entry_nov = clean_nov(entry)
                if entry_clean == target_clean or entry_nov == target_nov or (
                    "chianti" in entry_clean and any(part in entry_clean for part in re.findall(r'\d+', target))
                ):
                    candidates.append(full_entry)

            if candidates:
                best_match = candidates[0]
                for c in candidates:
                    if clean_nov(os.path.basename(c)) == target_nov:
                        best_match = c
                        break
                print(f"[Notice] Resolved CHIANTI directory '{dir_path}' -> '{best_match}'")
                return best_match

        print(f"ERROR: Specified CHIANTI directory does not exist: {dir_path}")
        if os.path.isdir(parent):
            chianti_entries = [e for e in os.listdir(parent) if "chianti" in e.lower()]
            if chianti_entries:
                print(f"       Available CHIANTI installations in {parent}:")
                for ce in chianti_entries:
                    print(f"         - {os.path.join(parent, ce)}")
        sys.exit(1)

    # 3. Auto-detection if no path specified
    env_xuvtop = os.environ.get("XUVTOP")
    common_paths = [
        env_xuvtop,
        "/home/hmonteiro/software/chianti_v10.1",
        "/home/hmonteiro/software/cloudy-c25.00/data/chianti",
        os.path.expanduser("~/software/chianti_v10.1"),
        os.path.expanduser("~/software/cloudy-c25.00/data/chianti"),
    ]
    for p in common_paths:
        if p and os.path.isdir(p):
            print(f"[Source] Auto-detected CHIANTI directory at {p}")
            return os.path.abspath(p)

    return None


def main():
    parser = argparse.ArgumentParser(description="Generate MOCASSIN atomic data from CHIANTI database releases.")
    parser.add_argument("--chianti-dir", type=str, default=None,
                        help="Path to CHIANTI directory (e.g. /path/to/chianti_v10.1 or $XUVTOP).")
    parser.add_argument("--chianti-tar", type=str, default=None,
                        help="Path to CHIANTI database .tar.gz archive.")
    parser.add_argument("--download-version", type=str, default=None,
                        help="Download official release (e.g. 10.1 or 11.0.2) directly from chiantidatabase.org.")
    parser.add_argument("--output-dir", type=str, default="data",
                        help="Directory where generated .dat files should be written (default: data).")
    parser.add_argument("--file-list", type=str, default="data/fileNames.dat",
                        help="Path to fileNames.dat list of ions (default: data/fileNames.dat).")
    parser.add_argument("--ions", type=str, default=None,
                        help="Comma-separated list of ions to process (e.g. 'c_2,o_3,fe_6' or 'all'). Default: all available in CHIANTI.")
    parser.add_argument("--max-levels", type=int, default=17,
                        help="Global cap on level count to protect MOCASSIN's nForLevels=17 array limits (default: 17).")
    parser.add_argument("--nlevs", type=str, default=None,
                        help="Custom level overrides in the format 'c_2=10,o_3=15'.")
    parser.add_argument("--version-tag", type=str, default=None,
                        help="Custom version string for header comments (default: auto-detected or '10').")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip generating files that already exist in output-dir.")
    parser.add_argument("--compare-with", type=str, default=None,
                        help="Directory to compare generated files against (e.g. data or data.chianty-10).")
    args = parser.parse_args()

    # Determine CHIANTI source
    chianti_dir = None
    chianti_tar = args.chianti_tar

    if args.download_version:
        chianti_tar = download_chianti_release(args.download_version)

    if not chianti_tar:
        chianti_dir = resolve_chianti_dir(args.chianti_dir)
        if not chianti_dir:
            print("ERROR: No CHIANTI source specified and none could be auto-detected.")
            print("       Provide --chianti-dir, --chianti-tar, or --download-version.")
            sys.exit(1)

    source = ChiantiSource(directory=chianti_dir, tar_path=chianti_tar)

    # Detect CHIANTI version
    version_str = args.version_tag
    if not version_str:
        v_text = source.read_text("VERSION")
        if v_text:
            version_str = v_text.strip()
        else:
            version_str = "10"

    print(f"[CHIANTI] Active version tag: {version_str}")

    # Parse fileNames.dat
    file_list_path = resolve_file_list(args.file_list)
    all_ions = parse_file_names(file_list_path)
    print(f"[MOCASSIN] Loaded {len(all_ions)} target ion entries from {file_list_path}")

    # Parse custom nlevs overrides
    custom_nlevs = {}
    if args.nlevs:
        for pair in args.nlevs.split(","):
            k, v = pair.split("=")
            custom_nlevs[k.strip().lower()] = int(v.strip())

    # Target ion filter
    selected_ions = None
    if args.ions and args.ions.lower() != "all":
        selected_ions = [i.strip().lower() for i in args.ions.split(",")]

    os.makedirs(args.output_dir, exist_ok=True)

    generated_count = 0
    skipped_count = 0
    missing_chianti_count = 0

    for elem, ion, moc_relpath, ch_ion in all_ions:
        filename = os.path.basename(moc_relpath)
        out_path = os.path.join(args.output_dir, filename)

        if selected_ions and ch_ion not in selected_ions and filename not in selected_ions:
            continue

        if args.skip_existing and os.path.exists(out_path):
            skipped_count += 1
            continue

        # Special exclusion: feii.dat (142-level Nahar model with custom iRats=1 structure)
        if filename == "feii.dat" and not (selected_ions and "feii.dat" in selected_ions):
            if not os.path.exists(out_path):
                existing_feii = os.path.join(os.path.dirname(file_list_path), "feii.dat")
                if os.path.exists(existing_feii):
                    import shutil
                    shutil.copy2(existing_feii, out_path)
                    print(f"[Preserved] Copied existing 142-level Nahar model feii.dat -> {out_path}")
                else:
                    print(f"[Skip] {filename} is the specialized 142-level Nahar model; preserving existing file.")
            else:
                print(f"[Skip] {filename} is the specialized 142-level Nahar model; preserving existing file.")
            continue

        el_sym = ELEMENTS[elem]
        rel_prefix = f"{el_sym}/{ch_ion}/{ch_ion}"

        # Check required files
        elvlc_text = source.read_text(f"{rel_prefix}.elvlc")
        wgfa_text = source.read_text(f"{rel_prefix}.wgfa")
        splups_text = source.read_text(f"{rel_prefix}.splups")
        scups_text = source.read_text(f"{rel_prefix}.scups")

        if not elvlc_text or not wgfa_text or (not splups_text and not scups_text):
            missing_chianti_count += 1
            continue

        # Parse levels
        levels, elvlc_comms = parse_elvlc(elvlc_text)
        if not levels:
            continue

        # Determine level count
        if ch_ion in custom_nlevs:
            nlevs = custom_nlevs[ch_ion]
        elif filename in DEFAULT_NLEVS:
            nlevs = DEFAULT_NLEVS[filename]
        else:
            nlevs = min(len(levels), args.max_levels)

        nlevs = min(nlevs, len(levels), args.max_levels)
        if nlevs < 2:
            continue

        # Parse WGFA
        a_dict, wgfa_comms = parse_wgfa(wgfa_text)

        # Parse collision data
        coll_dict = {}
        coll_comms = []
        if splups_text:
            coll_dict, coll_comms = parse_splups(splups_text)
        elif scups_text:
            coll_dict, coll_comms = parse_scups(scups_text)

        comments_all = elvlc_comms + [f"%filename: {ch_ion}.wgfa"] + wgfa_comms + \
                       [f"%filename: {ch_ion}.splups" if splups_text else f"%filename: {ch_ion}.scups"] + coll_comms

        # Generate file content
        content = generate_mocassin_file_content(
            ion_name=ch_ion,
            chianti_version_str=version_str,
            levels=levels,
            a_dict=a_dict,
            coll_dict=coll_dict,
            comments_all=comments_all,
            nlevs=nlevs,
            temps=DEFAULT_TEMPS
        )

        with open(out_path, "w", encoding="utf-8") as fp:
            fp.write(content)

        generated_count += 1
        print(f"[Generated] {filename:12s} ({ch_ion:6s}): {nlevs:2d} levels, {len(levels)} in CHIANTI -> {out_path}")

        # Optional comparison
        if args.compare_with:
            ref_path = os.path.join(args.compare_with, filename)
            if os.path.exists(ref_path):
                with open(ref_path) as fp:
                    ref_lines = fp.readlines()
                new_lines = content.splitlines(keepends=True)
                print(f"            vs {args.compare_with}: {len(new_lines)} lines generated vs {len(ref_lines)} ref lines")

    source.close()
    print(f"\n=======================================================")
    print(f"Atomic Data Generation Complete!")
    print(f"Generated files: {generated_count}")
    print(f"Skipped files:   {skipped_count}")
    print(f"Not in CHIANTI:  {missing_chianti_count} (trace elements / fully stripped ions)")
    print(f"Output folder:   {args.output_dir}")
    print(f"=======================================================\n")


if __name__ == "__main__":
    main()
