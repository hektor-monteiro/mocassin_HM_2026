#!/usr/bin/env python3
"""
test_hei_converter.py

Comprehensive verification suite for the MOCASSIN Porter He I converter:
1. Verifies dataset discovery and extraction from HDF5.
2. Tests 3-parameter non-linear curve fitting and checks residual thresholds.
3. Tests exact 102-line formatting and column structures.
4. Compiles and executes a standalone Fortran 90 test harness to verify
   runtime compatibility with source/hydro_mod.f90 (readHeIRecLines).
5. Verifies physical consistency against raw Porter tabulations.
"""

import os
import sys
import tempfile
import subprocess
import numpy as np

# Add parent directory to path
tools_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, tools_dir)

from generate_porter_hei_data import (
    locate_porter_hdf5,
    generate_porter_coefficients,
    write_mocassin_data_file,
    fit_emissivity_formula,
    MOCASSIN_WAVELENGTHS,
    REFERENCE_DENSITIES
)


def test_hdf5_discovery():
    print("[TEST 1/5] Testing Porter et al. HDF5 dataset discovery...")
    path = locate_porter_hdf5()
    assert path is not None, "Failed to locate he_i_rec_Pal12-Pal13.hdf5 on system."
    assert os.path.isfile(path), f"Discovered path does not exist: {path}"
    print(f"  -> Successfully located dataset at: {path}")


def test_curve_fitting():
    print("[TEST 2/5] Testing non-linear curve fitting residuals...")
    path = locate_porter_hdf5()
    results = generate_porter_coefficients(path, verbose=False)
    assert len(results) == 102, f"Expected 102 fits (34 lines x 3 densities), got {len(results)}"

    residuals = [v[4] for v in results.values()]
    median_res = np.median(residuals)
    max_res = np.max(residuals)

    print(f"  -> Median residual: {median_res:.2f}% (threshold: < 3.0%)")
    print(f"  -> Maximum residual: {max_res:.2f}% (threshold: < 15.0%)")
    assert median_res < 3.0, f"Median residual {median_res:.2f}% exceeds threshold 3.0%"
    assert max_res < 15.0, f"Maximum residual {max_res:.2f}% exceeds threshold 15.0%"
    print("  -> Curve fitting accuracy verified.")


def test_file_structure():
    print("[TEST 3/5] Testing output file structure (102 lines, 4 columns)...")
    path = locate_porter_hdf5()
    results = generate_porter_coefficients(path, verbose=False)

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        write_mocassin_data_file(results, tmp_path)
        with open(tmp_path) as f:
            lines = [l.strip() for l in f if l.strip()]

        assert len(lines) == 102, f"Expected exactly 102 lines, got {len(lines)}"

        for idx, line in enumerate(lines):
            parts = line.split()
            assert len(parts) == 4, f"Line {idx+1} does not have 4 columns: {line}"
            A, b, c, wl = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
            assert A > 0.0, f"Line {idx+1}: coefficient A must be positive, got {A}"
            assert wl > 2000.0 and wl < 25000.0, f"Line {idx+1}: unexpected wavelength {wl}"

        print("  -> File formatting and column structure verified.")
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_fortran_reader_compatibility():
    print("[TEST 4/5] Testing Fortran reader runtime compatibility...")
    path = locate_porter_hdf5()
    results = generate_porter_coefficients(path, verbose=False)

    with tempfile.TemporaryDirectory() as tmpdir:
        data_file = os.path.join(tmpdir, "HeIrecLines.dat")
        write_mocassin_data_file(results, data_file)

        fortran_source = f"""
program test_he1_runtime
  implicit none
  real :: HeIrecLineCoeff(34, 3, 4)
  integer :: iline, iden, j, ios
  character(len=256) :: fname

  fname = '{data_file}'
  open(unit=13, file=trim(fname), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'ERROR: Cannot open HeIrecLines.dat'
    stop 1
  end if

  do iden = 1, 3
     do iline = 1, 34
        read(13, *, iostat=ios) (HeIrecLineCoeff(iline, iden, j), j = 1, 4)
        if (ios /= 0) then
           print *, 'ERROR: Parsing error at iden=', iden, ' iline=', iline
           stop 2
        end if
        HeIrecLineCoeff(iline, iden, 1) = HeIrecLineCoeff(iline, iden, 1) * 1.e25
     end do
  end do
  close(13)

  ! Check that key lines have valid positive values
  if (HeIrecLineCoeff(16, 1, 1) <= 0.0 .or. HeIrecLineCoeff(16, 1, 4) < 5800.0) then
     print *, 'ERROR: He I 5876 corrupted'
     stop 3
  end if

  print *, 'FORTRAN_SUCCESS'
end program test_he1_runtime
"""
        f90_file = os.path.join(tmpdir, "test.f90")
        bin_file = os.path.join(tmpdir, "test_bin")

        with open(f90_file, "w") as f:
            f.write(fortran_source)

        compile_res = subprocess.run(["gfortran", "-O2", f90_file, "-o", bin_file],
                                     capture_output=True, text=True)
        assert compile_res.returncode == 0, f"Fortran compilation failed:\n{compile_res.stderr}"

        run_res = subprocess.run([bin_file], capture_output=True, text=True)
        assert run_res.returncode == 0, f"Fortran execution failed:\n{run_res.stderr}"
        assert "FORTRAN_SUCCESS" in run_res.stdout, "Fortran execution did not complete cleanly."
        print("  -> Fortran reader runtime compatibility verified (0 errors).")


def test_physical_consistency():
    print("[TEST 5/5] Testing physical consistency against raw Porter tabulations...")
    import h5py
    path = locate_porter_hdf5()
    results = generate_porter_coefficients(path, verbose=False)

    with h5py.File(path, "r") as f:
        data = f["updated_data"][:]

    # Test at Ne = 100 cm^-3 (log Ne = 2.0), Te = 10000 K (T4 = 1.0)
    sub = data[(data["DENS"] == 2.0) & (data["TEMP"] == 10000.0)][0]

    check_lines = [
        (15, "5876.0", 5875.66),
        (10, "4471.0", 4471.50),
        (16, "6678.0", 6678.16),
        (20, "10830.0", 10830.25),
    ]

    for iline, field, wl in check_lines:
        raw_val = sub[field]
        A, b, c, _, _ = results[(iline, 1)]
        fit_val = fit_emissivity_formula(1.0, A, b, c)
        diff_pct = abs(fit_val - raw_val) / raw_val * 100.0
        print(f"  -> {wl:8.2f} A at 10,000 K: raw={raw_val:.4e}, fit={fit_val:.4e}, diff={diff_pct:.2f}%")
        assert diff_pct < 2.0, f"Line {wl:.2f} A fit diff {diff_pct:.2f}% exceeds 2.0%"

    print("  -> Physical consistency across benchmark lines verified.")


def main():
    print("=" * 60)
    print("MOCASSIN Porter He I Converter Verification Suite")
    print("=" * 60)
    test_hdf5_discovery()
    test_curve_fitting()
    test_file_structure()
    test_fortran_reader_compatibility()
    test_physical_consistency()
    print("=" * 60)
    print("ALL TESTS PASSED SUCCESSFULLY!")
    print("=" * 60)


if __name__ == "__main__":
    main()
