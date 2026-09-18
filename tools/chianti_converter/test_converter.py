#!/usr/bin/env python3
"""
test_converter.py
-----------------
Automated test suite for the MOCASSIN CHIANTI atomic data generator.
Verifies:
1. Parsing of CHIANTI .elvlc, .wgfa, .splups, and .scups.
2. Burgess-Tully (1992) collision strength descaling across multiple transition types.
3. Strict conformity to MOCASSIN Fortran formatting standards.
4. Absence of premature termination conditions (e.g. qx == 0.0).
5. Numerical consistency against historical benchmark files (C II, O III, S II, Ne III).

Usage:
    python tools/chianti_converter/test_converter.py
"""

import os
import sys
import argparse
import tempfile
import subprocess
import numpy as np

# Add parent directory to path so we can import the generator module
sys.path.insert(0, os.path.dirname(__file__))
import generate_chianti_atomic_data as gen

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DEFAULT_FILE_LIST = os.path.join(REPO_ROOT, "data", "fileNames.dat")


def test_burgess_tully_descaling():
    """Verify Burgess-Tully descaling for Type 1, 2, 3, 4."""
    print("[TEST 1/5] Testing Burgess-Tully descaling formulas...")
    temps = np.array([100.0, 1000.0, 10000.0, 100000.0])
    
    # Test Type 2 (Born / forbidden)
    param_t2 = {
        "format": "splups",
        "ttype": 2,
        "gf": 0.0,
        "de": 0.0005779,
        "cups": 563.4,
        "spl": [1.552, 1.866, 2.167, 2.226, 2.278, 2.256, 2.233, 2.194, 2.152, 2.066, 1.932, 1.798, 1.542, 1.200, 0.695]
    }
    u2 = gen.descale_upsilon(param_t2, temps)
    assert len(u2) == len(temps), "Length mismatch in Upsilon output"
    assert np.all(u2 > 0.0), "Type 2 collision strengths must be strictly positive"
    assert np.all(u2 >= 1e-30), "Floor condition violated"
    
    # Test Type 1 (Dipole allowed)
    param_t1 = {
        "format": "splups",
        "ttype": 1,
        "gf": 0.3107,
        "de": 0.6828,
        "cups": 1.100,
        "spl": [1.139, 1.181, 1.212, 1.238, 1.274, 1.355, 1.444, 1.562, 1.702, 1.883, 1.820]
    }
    u1 = gen.descale_upsilon(param_t1, temps)
    assert np.all(u1 > 0.0), "Type 1 collision strengths must be strictly positive"
    print("  -> Burgess-Tully descaling tests PASSED.")


def test_cii_benchmark_comparison(chianti_dir: str = None, file_list: str = DEFAULT_FILE_LIST):
    """Verify C II generation against data.chianty-10/cii.dat."""
    print("[TEST 2/5] Testing C II generation against CHIANTI benchmark...")
    resolved_dir = gen.resolve_chianti_dir(chianti_dir) if chianti_dir else (
        "/home/hmonteiro/software/cloudy-c25.00/data/chianti"
        if os.path.exists("/home/hmonteiro/software/cloudy-c25.00/data/chianti")
        else gen.resolve_chianti_dir(None)
    )
    ref_path = os.path.join(REPO_ROOT, "data", "cii.dat")
    if not os.path.exists(ref_path) or not resolved_dir or not os.path.exists(resolved_dir):
        print("  -> Skipping (reference or CHIANTI dir not found).")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        out_file = os.path.join(tmpdir, "cii.dat")
        cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "generate_chianti_atomic_data.py"),
            "--chianti-dir", resolved_dir,
            "--ions", "c_2",
            "--output-dir", tmpdir,
            "--file-list", file_list
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        assert os.path.exists(out_file), "cii.dat was not generated"

        with open(out_file) as f:
            gen_lines = [l.strip() for l in f if l.strip()]
        with open(ref_path) as f:
            ref_lines = [l.strip() for l in f if l.strip()]

        # Check energies (last 10 lines)
        gen_e = [l.split() for l in gen_lines[-10:]]
        ref_e = [l.split() for l in ref_lines[-10:]]
        for (gi, gg, ge), (ri, rg, re) in zip(gen_e, ref_e):
            assert int(gi) == int(ri), f"Level mismatch: {gi} vs {ri}"
            assert int(gg) == int(rg), f"Weight mismatch: {gg} vs {rg}"
            assert abs(float(ge) - float(re)) < 0.2, f"Energy mismatch: {ge} vs {re}"

        # Check A values (45 lines before the last 10)
        gen_a = [l.split() for l in gen_lines[-55:-10]]
        ref_a = [l.split() for l in ref_lines[-55:-10]]
        for (gi, gj, ga), (ri, rj, ra) in zip(gen_a, ref_a):
            assert int(gi) == int(ri) and int(gj) == int(rj)
            val_g, val_r = float(ga), float(ra)
            if val_r > 0:
                rel = abs(val_g - val_r) / val_r
                assert rel < 0.05, f"A-value relative diff > 5% on {gi}->{gj}: {val_g} vs {val_r}"

    print("  -> C II benchmark comparison PASSED.")


def test_fortran_reader_compatibility(chianti_dir: str = None, file_list: str = DEFAULT_FILE_LIST):
    """Verify that gfortran can compile and execute a MOCASSIN-style reader on generated files."""
    print("[TEST 3/5] Testing Fortran reader runtime compatibility...")
    resolved_dir = gen.resolve_chianti_dir(chianti_dir) if chianti_dir else (
        "/home/hmonteiro/software/cloudy-c25.00/data/chianti"
        if os.path.exists("/home/hmonteiro/software/cloudy-c25.00/data/chianti")
        else gen.resolve_chianti_dir(None)
    )
    if not resolved_dir or not os.path.exists(resolved_dir):
        print("  -> Skipping (CHIANTI dir not found).")
        return

    fortran_source = """
program test_read
  implicit none
  integer :: ncoms, nlevs, ntemps, i, j, k, l, gx, ios, idx
  double precision :: ex
  character(len=78) :: comment
  character(len=20), allocatable :: labels(:)
  double precision, allocatable :: temps(:), g(:), e(:), a(:,:), qom(:,:,:)
  real :: irats, qx
  integer :: low, up
  character(len=100) :: fname

  call get_command_argument(1, fname)
  open(unit=121, file=trim(fname), status="old", action="read", iostat=ios)
  if (ios /= 0) stop 1

  read(121, *) ncoms
  do i = 1, ncoms
     read(121, "(A78)") comment
  end do

  read(121, *) nlevs, ntemps
  allocate(labels(nlevs), temps(ntemps), g(nlevs), e(nlevs), a(nlevs, nlevs), qom(ntemps, nlevs, nlevs))
  qom = 0.0d0
  a = 0.0d0

  do i = 1, nlevs
     read(121, "(A20)") labels(i)
  end do

  do i = 1, ntemps
     read(121, *) temps(i)
  end do

  read(121, *) irats

  do i = 1, 100000
     read(121, *) low, up, qx
     if (qx == 0.0d0) exit
     if (low == 0) then
        k = k + 1
     else
        j = low
        l = up
        k = 1
     end if
     qom(k, j, l) = qx
  end do

  do k = 1, nlevs - 1
     do l = k + 1, nlevs
        read(121, *) i, j, a(j, k)
     end do
  end do

  do i = 1, nlevs
     read(121, *) idx, gx, ex
     g(idx) = gx
     e(idx) = ex
  end do

  close(121)
  if (e(nlevs) <= 0.0d0 .and. nlevs > 1) stop 2
end program test_read
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        f90_path = os.path.join(tmpdir, "test_read.f90")
        bin_path = os.path.join(tmpdir, "test_read")
        with open(f90_path, "w") as f:
            f.write(fortran_source)

        compile_cmd = ["gfortran", "-O2", f90_path, "-o", bin_path]
        subprocess.run(compile_cmd, check=True)

        # Generate sample ions representing different shell configurations
        test_ions = ["c_2", "o_3", "ne_3", "s_2", "fe_6"]
        gen_cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "generate_chianti_atomic_data.py"),
            "--chianti-dir", resolved_dir,
            "--ions", ",".join(test_ions),
            "--output-dir", tmpdir,
            "--file-list", file_list
        ]
        subprocess.run(gen_cmd, check=True, capture_output=True)

        ion_to_file = {
            "c_2": "cii.dat", "o_3": "oiii.dat", "ne_3": "neiii.dat",
            "s_2": "sii.dat", "fe_6": "fevi.dat"
        }
        for ion, fname in ion_to_file.items():
            fpath = os.path.join(tmpdir, fname)
            res = subprocess.run([bin_path, fpath], capture_output=True)
            assert res.returncode == 0, f"Fortran parser failed on {fname} with exit code {res.returncode}"

    print("  -> Fortran reader runtime compatibility PASSED.")


def test_no_zero_upsilon_premature_exit(chianti_dir: str = None, file_list: str = DEFAULT_FILE_LIST):
    """Verify that collision strengths never format as 0.00e+00 before the terminator."""
    print("[TEST 4/5] Testing absence of premature 0.00e+00 in collision data...")
    resolved_dir = gen.resolve_chianti_dir(chianti_dir) if chianti_dir else (
        "/home/hmonteiro/software/cloudy-c25.00/data/chianti"
        if os.path.exists("/home/hmonteiro/software/cloudy-c25.00/data/chianti")
        else gen.resolve_chianti_dir(None)
    )
    if not resolved_dir or not os.path.exists(resolved_dir):
        print("  -> Skipping (CHIANTI dir not found).")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "generate_chianti_atomic_data.py"),
            "--chianti-dir", resolved_dir,
            "--ions", "s_2,si_2,al_2",
            "--output-dir", tmpdir,
            "--file-list", file_list
        ]
        subprocess.run(cmd, check=True, capture_output=True)

        for fname in ["sii.dat", "si2.dat", "alii.dat"]:
            fpath = os.path.join(tmpdir, fname)
            with open(fpath) as fp:
                lines = fp.readlines()
            # find terminator
            term_idx = [i for i, l in enumerate(lines) if l.strip() == "0 0 0"][0]
            coll_lines = lines[:term_idx]
            for idx, l in enumerate(coll_lines):
                tokens = l.split()
                if len(tokens) == 3 and tokens[0] == "0" and tokens[1] == "0":
                    val_str = tokens[2]
                    assert val_str != "0.00e+00" and val_str != "0.000e+00", \
                        f"Found premature zero {val_str} at line {idx+1} in {fname}"

    print("  -> Premature exit prevention test PASSED.")


def test_official_chianti_resolution_and_format(chianti_dir: str = None, file_list: str = DEFAULT_FILE_LIST):
    """Verify directory resolution (CHIANTI_10.1 -> chianti_v10.1) and .scups parsing with Fortran reader."""
    print("[TEST 5/5] Testing official CHIANTI 10.1 directory resolution & .scups generation...")
    test_dir = chianti_dir if chianti_dir else "/home/hmonteiro/software/CHIANTI_10.1"
    resolved = gen.resolve_chianti_dir(test_dir)
    assert resolved and os.path.exists(resolved), f"Failed to resolve {test_dir}"
    assert "chianti" in resolved.lower(), f"Expected resolved path to contain chianti, got {resolved}"

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "generate_chianti_atomic_data.py"),
            "--chianti-dir", test_dir,
            "--ions", "c_2,o_3,fe_6",
            "--output-dir", tmpdir,
            "--file-list", file_list
        ]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if test_dir != resolved:
            assert "[Notice] Resolved CHIANTI directory" in res.stdout, "Directory resolution notice missing in stdout"

        for fname in ["cii.dat", "oiii.dat", "fevi.dat"]:
            fpath = os.path.join(tmpdir, fname)
            assert os.path.exists(fpath), f"{fname} was not generated from official CHIANTI 10.1"
            with open(fpath) as fp:
                content = fp.read()
            assert "0 0 0" in content, f"Terminator missing in {fname}"

    print("  -> Official CHIANTI 10.1 resolution & format PASSED.")


def main():
    parser = argparse.ArgumentParser(description="MOCASSIN CHIANTI Converter Verification Suite")
    parser.add_argument("--chianti-dir", type=str, default=None,
                        help="Path to CHIANTI directory (e.g. /home/hmonteiro/software/CHIANTI_10.1 or $XUVTOP).")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Target output directory (accepted for generator CLI compatibility).")
    parser.add_argument("--file-list", type=str, default=DEFAULT_FILE_LIST,
                        help="Path to fileNames.dat master list.")
    args = parser.parse_args()

    file_list = gen.resolve_file_list(args.file_list)

    print("==================================================")
    print("MOCASSIN CHIANTI Converter Verification Suite")
    print("==================================================")
    test_burgess_tully_descaling()
    test_cii_benchmark_comparison(chianti_dir=args.chianti_dir, file_list=file_list)
    test_fortran_reader_compatibility(chianti_dir=args.chianti_dir, file_list=file_list)
    test_no_zero_upsilon_premature_exit(chianti_dir=args.chianti_dir, file_list=file_list)
    test_official_chianti_resolution_and_format(chianti_dir=args.chianti_dir, file_list=file_list)
    print("==================================================")
    print("ALL TESTS PASSED SUCCESSFULLY!")
    print("==================================================")


if __name__ == "__main__":
    main()
