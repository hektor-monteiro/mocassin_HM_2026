#!/usr/bin/env python3
"""
generate_porter_hei_data.py

Generates modern He I recombination line coefficients for MOCASSIN from the
Porter et al. (2012, MNRAS, 425, L28; erratum 2013, MNRAS, 433, L89) dataset.

Fits the standard 3-parameter formula:
    4*pi*j / (Ne * N(He+)) = A * (T4**b) * exp(c / T4)  [10^-25 erg*cm^3/s]
where T4 = Te / 1e4 K, across 3 reference densities (Ne = 10^2, 10^4, 10^6 cm^-3)
for all 34 MOCASSIN optical and near-infrared transitions.

The output format is 100% compatible with MOCASSIN's Fortran reader in
source/hydro_mod.f90 (subroutine readHeIRecLines).
"""

import os
import sys
import argparse
import numpy as np
from scipy.optimize import curve_fit

# Standard 34 transitions in MOCASSIN order
MOCASSIN_WAVELENGTHS = [
    4471.50, 2945.10, 3188.74, 3613.64, 3888.65, 3964.73, 4026.21,
    4120.82, 4387.93, 4437.55, 4471.50, 4713.17, 4921.93, 5015.68,
    5047.74, 5875.66, 6678.16, 7065.25, 7281.35, 9463.58, 10830.25,
    11969.06, 12527.49, 12784.92, 12790.50, 12968.43, 15083.65,
    17002.40, 18685.33, 18697.21, 19089.36, 19543.19, 20581.28, 21120.12
]

# Density reference points (Ne in cm^-3, log10(Ne))
REFERENCE_DENSITIES = [
    (1, 2.0),  # iden = 1: Ne = 10^2 cm^-3
    (2, 4.0),  # iden = 2: Ne = 10^4 cm^-3
    (3, 6.0),  # iden = 3: Ne = 10^6 cm^-3
]


def fit_emissivity_formula(T4, A, b, c):
    """Benjamin / MOCASSIN parametric emissivity formula."""
    return A * (T4 ** b) * np.exp(c / T4)


def locate_porter_hdf5(custom_path=None):
    """Locate the Porter et al. 2012/2013 HDF5 file."""
    candidates = []
    if custom_path:
        candidates.append(custom_path)

    # Check PyNeb package directory
    try:
        import pyneb
        pyneb_dir = os.path.dirname(pyneb.__file__)
        candidates.append(os.path.join(pyneb_dir, "atomic_data_fits", "he_i_rec_Pal12-Pal13.hdf5"))
    except ImportError:
        pass

    # Standard known paths in miniforge / conda environments
    candidates.extend([
        "/home/hmonteiro/miniforge3/lib/python3.13/site-packages/pyneb/atomic_data_fits/he_i_rec_Pal12-Pal13.hdf5",
        "/home/hmonteiro/miniforge3/lib/python3.12/site-packages/pyneb/atomic_data_fits/he_i_rec_Pal12-Pal13.hdf5",
        os.path.expanduser("~/.local/share/pyneb/atomic_data_fits/he_i_rec_Pal12-Pal13.hdf5"),
    ])

    for path in candidates:
        if os.path.isfile(path):
            return path

    return None


def generate_porter_coefficients(hdf5_path, verbose=False):
    """
    Extracts tabular data from the Porter HDF5 file and performs non-linear
    curve fits for each line and density.
    
    Returns:
        dict: (iline, iden) -> (A, b, c, wl, max_err_pct)
    """
    try:
        import h5py
    except ImportError:
        raise RuntimeError("h5py package is required to read the Porter HDF5 dataset.")

    if not os.path.isfile(hdf5_path):
        raise FileNotFoundError(f"Porter HDF5 file not found: {hdf5_path}")

    with h5py.File(hdf5_path, "r") as f:
        data = f["updated_data"][:]

    temps = np.unique(data["TEMP"])
    T4 = temps / 1e4

    p_fields = [name for name in data.dtype.names if name not in ["TEMP", "DENS"]]
    p_wls = [float(name) for name in p_fields]

    results = {}

    for iden_idx, log_ne in REFERENCE_DENSITIES:
        sub_data = data[data["DENS"] == log_ne]
        sort_idx = np.argsort(sub_data["TEMP"])
        sub_data = sub_data[sort_idx]

        for iline, mw in enumerate(MOCASSIN_WAVELENGTHS):
            # Find closest matching transition in Porter dataset
            pw = min(p_wls, key=lambda x: abs(x - mw))
            field_name = None
            for nm in p_fields:
                if abs(float(nm) - mw) < 5.0:
                    field_name = nm
                    break

            if field_name is None:
                raise ValueError(f"Could not find matching Porter transition for {mw:.2f} A")

            y = sub_data[field_name]

            # Fit 3 parameters: A, b, c
            p0 = [y[5], -1.0, 0.0]
            try:
                popt, _ = curve_fit(fit_emissivity_formula, T4, y, p0=p0, maxfev=10000)
            except Exception:
                # Fallback initial guess
                p0 = [y[0], -0.5, 0.1]
                popt, _ = curve_fit(fit_emissivity_formula, T4, y, p0=p0, maxfev=20000)

            A, b, c = popt[0], popt[1], popt[2]
            y_pred = fit_emissivity_formula(T4, A, b, c)
            max_err = np.max(np.abs((y_pred - y) / y)) * 100.0

            results[(iline, iden_idx)] = (A, b, c, mw, max_err)

            if verbose and (iline in [0, 6, 15, 16, 20]):
                print(f"  [iden={iden_idx}] {mw:8.2f} A: A={A:10.3e}, b={b:6.3f}, c={c:6.3f}, max_err={max_err:4.1f}%")

    return results


def write_mocassin_data_file(results, output_path):
    """
    Writes the exact 102-line Fortran data file:
        do iden = 1, 3
           do iline = 1, 34
              read(13,*) (HeIrecLineCoeff(iline,iden,j), j = 1, 4)
    """
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with open(output_path, "w") as f:
        for iden_idx, _ in REFERENCE_DENSITIES:
            for iline in range(len(MOCASSIN_WAVELENGTHS)):
                A, b, c, mw, _ = results[(iline, iden_idx)]
                f.write(f"{A:9.2e}\t{b:6.3f}\t{c:6.3f}\t {mw:7.2f}\n")

    print(f"[Success] Written 102 lines to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate modern He I recombination line coefficients for MOCASSIN (Porter et al. 2012/2013)."
    )
    parser.add_argument("-o", "--output", default="data/HeIrecLines.dat",
                        help="Target output file path (default: data/HeIrecLines.dat)")
    parser.add_argument("-s", "--source", default=None,
                        help="Path to Porter et al. HDF5 dataset (auto-detected if omitted)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print detailed fit diagnostics per line")
    args = parser.parse_args()

    hdf5_path = locate_porter_hdf5(args.source)
    if not hdf5_path:
        print("Error: Could not locate Porter et al. dataset ('he_i_rec_Pal12-Pal13.hdf5').", file=sys.stderr)
        print("Please install PyNeb ('pip install pyneb') or specify the path with --source.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading Porter et al. (2012/2013) dataset from: {hdf5_path}")
    results = generate_porter_coefficients(hdf5_path, verbose=args.verbose)

    # Print summary statistics
    all_errs = [v[4] for v in results.values()]
    print(f"Fitting complete across 34 transitions and 3 densities:")
    print(f"  Median max fit residual: {np.median(all_errs):.2f}%")
    print(f"  90th percentile residual: {np.percentile(all_errs, 90):.2f}%")
    print(f"  Maximum residual:        {np.max(all_errs):.2f}%")

    write_mocassin_data_file(results, args.output)


if __name__ == "__main__":
    main()
