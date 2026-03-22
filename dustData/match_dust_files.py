#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 13:25:29 2025

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath # For complex mathematical operations like log

plt.close('all')

############################################################
# def interpolate_dust_properties(base_particle, target_particle):
#     """
#     Interpolates the optical constants (n and k) of a target dust particle
#     onto the wavelength grid of a base particle.

#     If the target particle's wavelength range is smaller than the base particle's,
#     the base particle's n and k values are used for the missing ranges.

#     Args:
#         base_particle (np.ndarray): A NumPy array with shape (N, 3)
#                                     representing [lambda, n, k] for the base particle.
#         target_particle (np.ndarray): A NumPy array with shape (M, 3)
#                                       representing [lambda, n, k] for the target particle to be interpolated.

#     Returns:
#         np.ndarray: A new NumPy array with shape (N, 3) containing the interpolated
#                     n and k values on the base particle's lambda grid.
#     """
#     # Unpack the arrays for clarity
#     lambda_base = base_particle[:, 0]
#     n_base = base_particle[:, 1]
#     k_base = base_particle[:, 2]

#     lambda_target = target_particle[:, 0]
#     n_target = target_particle[:, 1]
#     k_target = target_particle[:, 2]

#     # Perform linear interpolation of n and k values from the target particle
#     # onto the base particle's wavelength grid.
#     # For wavelengths outside the target's range, np.interp will use the
#     # specified 'left' and 'right' arguments. We set these to be the
#     # corresponding values from the base particle itself.
#     n_interpolated = np.interp(
#         lambda_base,
#         lambda_target,
#         n_target,
#         left=n_target[0],
#         right=n_target[-1]
#     )

#     k_interpolated = np.interp(
#         lambda_base,
#         lambda_target,
#         k_target,
#         left=k_target[0],
#         right=k_target[-1]
#     )
    
#     # check intervals and adjust accordingly
#     # if target min > base min
#     if (lambda_target[0] > lambda_base[0]):
#         ind = np.where(lambda_base < lambda_target[0])
#         n_interpolated[ind] = (n_base+n_interpolated)[ind]/2
#         k_interpolated[ind] = (k_base+k_interpolated)[ind]/2
        
#     # if target max < base max
#     if (lambda_target[-1] < lambda_base[-1]):
#         ind = np.where(lambda_base > lambda_target[-1])
#         n_interpolated[ind] = (n_base+n_interpolated)[ind]/2
#         k_interpolated[ind] = (k_base+k_interpolated)[ind]/2

#     # Combine the base lambda values with the new interpolated n and k values.
#     # The result is stacked horizontally to create the final [lambda, n, k] array.
#     n_interpolated = (n_interpolated)
#     interpolated_properties = np.vstack((lambda_base, n_interpolated, k_interpolated)).T

#     return interpolated_properties
###############################################################################
# def interpolate_dust_properties(base_particle: np.ndarray, target_particle: np.ndarray) -> np.ndarray:
#     """
#     Interpolates the optical constants (n and k) of a target dust particle
#     onto the wavelength grid of a base particle.

#     For wavelengths outside the target particle's range, the base particle's
#     original n and k values are used.

#     Args:
#         base_particle (np.ndarray): An array with shape (N, 3)
#                                     representing [lambda, n, k] for the base particle.
#         target_particle (np.ndarray): An array with shape (M, 3)
#                                       representing [lambda, n, k] for the target particle.

#     Returns:
#         np.ndarray: A new array with shape (N, 3) containing the interpolated
#                     n and k values on the base particle's lambda grid.
#     """
#     # Unpack the arrays for clarity
#     lambda_base = base_particle[:, 0]
#     n_base = base_particle[:, 1]
#     k_base = base_particle[:, 2]

#     lambda_target = target_particle[:, 0]
#     n_target = target_particle[:, 1]
#     k_target = target_particle[:, 2]

#     # Perform linear interpolation of n and k from the target onto the base grid
#     # For points outside the target's range, np.interp will use the endpoint values
#     n_interpolated = np.interp(lambda_base, lambda_target, n_target)
#     k_interpolated = np.interp(lambda_base, lambda_target, k_target)

#     # Identify the wavelengths in the base grid that are outside the target's range
#     extrapolation_mask = (lambda_base < lambda_target[0]) | (lambda_base > lambda_target[-1])

#     # In these regions, replace the interpolated values with the original base values
#     n_interpolated[extrapolation_mask] = n_base[extrapolation_mask]
#     k_interpolated[extrapolation_mask] = k_base[extrapolation_mask]

#     # Recombine the columns into the final [lambda, n, k] array
#     interpolated_properties = np.vstack((lambda_base, n_interpolated, k_interpolated)).T

#     return interpolated_properties

###############################################################################
def interpolate_dust_properties(base_particle, target_particle):
    """
    Interpolates optical constants (n, k) of a target particle onto a base grid
    with a smooth transition in non-overlapping regions.

    In regions where the target particle has no data, this function creates a
    smooth blend from the base particle's values to the target particle's
    endpoint values over a specified number of data points.

    Args:
        base_particle (np.ndarray): An array with shape (N, 3) for the base
                                    particle [lambda, n, k].
                                    Must be sorted by wavelength.
        target_particle (np.ndarray): An array with shape (M, 3) for the target
                                      particle [lambda, n, k].
                                      Must be sorted by wavelength.
        transition_points (int): The number of data points over which to blend
                                 the two datasets at the boundaries.

    Returns:
        np.ndarray: A new array with shape (N, 3) with interpolated and blended
                    n and k values on the base particle's lambda grid.
    """
    
    transition_points: int = 500
    
    # Unpack the arrays
    lambda_base = base_particle[:, 0]
    n_base = base_particle[:, 1]
    k_base = base_particle[:, 2]

    lambda_target = target_particle[:, 0]
    n_target = target_particle[:, 1]
    k_target = target_particle[:, 2]

    # Perform initial interpolation. This gives us the core interpolated values
    # and the constant-extrapolated values for the non-overlapping regions.
    n_interp = np.interp(lambda_base, lambda_target, n_target)
    k_interp = np.interp(lambda_base, lambda_target, k_target)

    # Create the final arrays, starting with the base particle's values.
    n_final = np.copy(n_base)
    k_final = np.copy(k_base)

    # Find the indices corresponding to the overlapping region
    # searchsorted finds where elements should be inserted to maintain order
    start_overlap_idx = np.searchsorted(lambda_base, lambda_target[0], side='left')
    end_overlap_idx = np.searchsorted(lambda_base, lambda_target[-1], side='right')

    # In the fully overlapping region, use the interpolated target values
    n_final[start_overlap_idx:end_overlap_idx] = n_interp[start_overlap_idx:end_overlap_idx]
    k_final[start_overlap_idx:end_overlap_idx] = k_interp[start_overlap_idx:end_overlap_idx]

    # --- Handle the smooth transition for shorter wavelengths (left side) ---
    if start_overlap_idx > 0 and transition_points > 0:
        # Define the transition zone
        trans_len = min(start_overlap_idx, transition_points)
        trans_start_idx = start_overlap_idx - trans_len
        zone = slice(trans_start_idx, start_overlap_idx)

        # Create blending weights (alpha) from 0 (pure base) to 1 (pure target)
        alpha = np.linspace(0, 1, trans_len)

        # Apply the weighted average
        n_final[zone] = (1 - alpha) * n_base[zone] + alpha * n_interp[zone]
        k_final[zone] = (1 - alpha) * k_base[zone] + alpha * k_interp[zone]

    # --- Handle the smooth transition for longer wavelengths (right side) ---
    if end_overlap_idx < len(lambda_base) and transition_points > 0:
        # Define the transition zone
        trans_len = min(len(lambda_base) - end_overlap_idx, transition_points)
        trans_end_idx = end_overlap_idx + trans_len
        zone = slice(end_overlap_idx, trans_end_idx)

        # Create blending weights (alpha) from 1 (pure target) to 0 (pure base)
        alpha = np.linspace(1, 0, trans_len)

        # Apply the weighted average
        n_final[zone] = (1 - alpha) * n_base[zone] + alpha * n_interp[zone]
        k_final[zone] = (1 - alpha) * k_base[zone] + alpha * k_interp[zone]

    # Recombine the columns into the final array
    interpolated_properties = np.vstack((lambda_base, n_final, k_final)).T

    return interpolated_properties
###############################################################################
def get_effective_optical_constants(n_a, k_a, wa, n_b, k_b, wb, n_c, k_c, wc):
    """
    Calculates the effective isotropic optical constants (n_eff, k_eff) for a
    material with three distinct polarization axes, assuming random orientation.

    This is done by averaging the complex dielectric function. This method is
    standard for modeling anisotropic dust grains in astrophysics.

    Args:
        n_a (np.ndarray): Array of refractive indices for the 'a' axis.
        k_a (np.ndarray): Array of extinction coefficients for the 'a' axis.
        wa  (np.ndarray): Array of weights for the 'a' axis.
        n_b (np.ndarray): Array of refractive indices for the 'b' axis.
        k_b (np.ndarray): Array of extinction coefficients for the 'b' axis.
        wb  (np.ndarray): Array of weights for the 'b' axis.
        n_c (np.ndarray): Array of refractive indices for the 'c' axis.
        k_c (np.ndarray): Array of extinction coefficients for the 'c' axis.
        wc  (np.ndarray): Array of weights for the 'c' axis.
        All input arrays must have the same length (i.e., correspond to the
        same wavelength points).

    Returns:
        tuple[np.ndarray, np.ndarray]: A tuple containing two numpy arrays:
            - n_eff: The effective refractive index.
            - k_eff: The effective extinction coefficient.
    """
    # --- Step 1: Convert n and k to the complex dielectric function (epsilon) ---
    # The dielectric function is epsilon = epsilon_1 + i*epsilon_2
    # where epsilon_1 = n^2 - k^2 and epsilon_2 = 2nk.
    # We can represent this directly using Python's complex numbers.
    
    epsilon_a = (n_a**2 - k_a**2) + 1j * (2 * n_a * k_a)
    epsilon_b = (n_b**2 - k_b**2) + 1j * (2 * n_b * k_b)
    epsilon_c = (n_c**2 - k_c**2) + 1j * (2 * n_c * k_c)

    # --- Step 2: Average the dielectric functions for random orientation ---
    # For a randomly oriented collection of grains, the effective dielectric
    # function is the arithmetic mean of the functions along the three axes.
    
    epsilon_eff = wa*epsilon_a + wb*epsilon_b + wc*epsilon_c

    # Extract the real (eps_1) and imaginary (eps_2) parts of the effective function
    eps_eff_1 = np.real(epsilon_eff)
    eps_eff_2 = np.imag(epsilon_eff)

    # --- Step 3: Convert the effective dielectric function back to n and k ---
    # We use the inverse relations to find n_eff and k_eff.
    
    # Calculate the modulus of the effective dielectric function squared
    mod_epsilon_sq = eps_eff_1**2 + eps_eff_2**2
    
    n_eff = np.sqrt((np.sqrt(mod_epsilon_sq) + eps_eff_1) / 2.0)
    k_eff = np.sqrt((np.sqrt(mod_epsilon_sq) - eps_eff_1) / 2.0)

    # --- Step 4: Return the final effective optical constants ---
    return n_eff, k_eff
    

###############################################################################
# calculate coeficients
def calculate_mie_cabs(k, V, alpha):
    """
    Calculates Mie absorption cross-section (Cabs).
    Cabs = k * V * imag(3*alpha)
    """
    return k * V * alpha.imag * 3

def calculate_mie_csca(k, V, alpha):
    """
    Calculates Mie scattering cross-section (Csca).
    Csca = k^4 * V^2 * 3 * absolute(alpha)^2 / (2*pi)
    Note: Original email had 2*pi in denominator.
    """
    return (k**4) * (V**2) * 3 * (abs(alpha)**2) / (2 * np.pi)

def calculate_cde_cabs(k, V, alpha):
    """
    Calculates CDE absorption cross-section (Cabs).
    Cabs = k * V * alpha.imag
    """
    return k * V * alpha.imag

def calculate_cde_csca(k, V, alpha, sigma):
    """
    Calculates CDE scattering cross-section (Csca).
    Csca = k * V * imag(alpha) / sigma
    Handles potential division by zero, NaN, or Inf in sigma using an epsilon check.
    """
    # Define a small epsilon to check if sigma is effectively zero.
    # This helps catch floating-point numbers that are very close to zero but not exactly zero.
    EPSILON = 1e-18

    if np.isnan(sigma) or np.isinf(sigma) or abs(sigma) < EPSILON:
        # Return NaN if sigma is problematic (NaN, Inf, or effectively zero)
        # This will prevent the ZeroDivisionError and the NaN values will not be plotted.
        return np.nan
    return k * V * alpha.imag / sigma

def process_nk(wave_micron, n_values, k_values, a_radius_microns):
    # Convert particle radius to cm for volume calculation
    particle_radius_cm = a_radius_microns * 1e-4 # 1 micron = 1e-4 cm
    # Assuming spherical particle volume
    V = (4/3) * np.pi * (particle_radius_cm**3)
    # Calculate geometrical cross-section (Ag)
    Ag = np.pi * (particle_radius_cm**2)


    # Lists to store calculated cross-sections and efficiencies
    mie_cabs_list = []
    mie_csca_list = []
    cde_cabs_list = []
    cde_csca_list = []
    mie_qabs_list = [] # New list for Mie Absorption Efficiency
    cde_qabs_list = [] # New list for CDE Absorption Efficiency


    for i in range(len(wave_micron)):
        lam_micron = wave_micron[i]
        n_val = n_values[i]
        k_val = k_values[i]

        # Convert wavelength from microns to cm for wavenumber calculation
        lam_cm = lam_micron * 1e-4 # 1 micron = 1e-4 cm
        k_wavenumber = (2 * np.pi) / lam_cm # k is wavenumber (2pi/lambda(cm))

        # Complex refractive index 'm'
        m = complex(n_val, k_val)

        # --- Mie Theory Calculations ---
        # alpha for Mie theory
        mie_alpha = (m**2 - 1) / (m**2 + 2)
        mie_cabs = calculate_mie_cabs(k_wavenumber, V, mie_alpha)
        mie_csca = calculate_mie_csca(k_wavenumber, V, mie_alpha)
        mie_cabs_list.append(mie_cabs)
        mie_csca_list.append(mie_csca)
        # Calculate Mie Absorption Efficiency (Qabs)
        mie_qabs_list.append(mie_cabs / Ag if Ag != 0 else np.nan)


        # --- CDE Theory Calculations ---
        # alpha for CDE theory
        # Handle log(m^2) carefully, especially if m^2 is zero or negative real.
        # cmath.log handles principal value (n=0) by default.
        try:
            cde_alpha = (2 * m**2 * cmath.log(m**2)) / (m**2 - 1) - 2
        except ZeroDivisionError:
            print(f"Warning: Division by zero when calculating CDE alpha for m^2-1 at wavelength {lam_micron} micron. Setting alpha to NaN.")
            cde_alpha = complex(np.nan, np.nan)
        except Exception as e:
            print(f"Warning: Error calculating CDE alpha for wavelength {lam_micron} micron: {e}. Setting alpha to NaN.")
            cde_alpha = complex(np.nan, np.nan)

        # sigma for CDE theory
        # Handle potential zero in denominator abs(m^2 - 1)^2
        try:
            sigma_denominator = (abs(m**2 - 1))**2
            if sigma_denominator == 0:
                print(f"Warning: Division by zero when calculating CDE sigma for abs(m^2-1)^2 at wavelength {lam_micron} micron. Setting sigma to NaN.")
                sigma = np.nan
            else:
                sigma = (6 * np.pi / (V * k_wavenumber**3)) * (m**2).imag / sigma_denominator
        except ZeroDivisionError: # This might not catch all cases, explicit check added above
            print(f"Warning: Division by zero when calculating CDE sigma (V*k^3) at wavelength {lam_micron} micron. Setting sigma to NaN.")
            sigma = np.nan
        except Exception as e:
            print(f"Warning: Error calculating CDE sigma for wavelength {lam_micron} micron: {e}. Setting sigma to NaN.")
            sigma = np.nan

        cde_cabs = calculate_cde_cabs(k_wavenumber, V, cde_alpha)
        cde_csca = calculate_cde_csca(k_wavenumber, V, cde_alpha, sigma)
        cde_cabs_list.append(cde_cabs)
        cde_csca_list.append(cde_csca)
        # Calculate CDE Absorption Efficiency (Qabs)
        cde_qabs_list.append(cde_cabs / Ag if Ag != 0 else np.nan)


    return mie_cabs_list, mie_csca_list, cde_cabs_list, cde_csca_list, mie_qabs_list, cde_qabs_list, V, k_wavenumber # Return k_wavenumber for reference
    

###############################################################################
# Base file of Astrophysical Silicates

astrosilfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/astroSilD2003.nk'
astrosildata = np.loadtxt(astrosilfile,max_rows=837,skiprows=2)

plt.figure()
plt.xscale('log')
plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

# increase sampling of astro silicate data
new_lambda = np.logspace(np.log10(astrosildata[:,0].min()),np.log10(astrosildata[:,0].max()), 2000)
new_Sil_n = np.interp(new_lambda,astrosildata[:,0],astrosildata[:,1])
new_Sil_k = np.interp(new_lambda,astrosildata[:,0],astrosildata[:,2])


astrosildata = np.vstack((new_lambda,new_Sil_n,new_Sil_k)).T

plt.plot(astrosildata[:,0], astrosildata[:,1],'C1',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C1',label='astrosil k')

header = "nk\nSsil_dl_base  1400. 3.6  0.588 20.077"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/astroSilD2003_base.nk'
np.savetxt(output_file, astrosildata, header=header, comments='', fmt='%.6E')

###############################################################################
# Water ice from Warren and Brandt 2008: Ice; n,k 0.00443–2e6 µm
# S. G. Warren and R. E. Brandt. Optical constants of ice from the ultraviolet to the microwave: A revised compilation. J. Geophys. Res. 113, D14220 (2008)

h2ofile08 = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_Warren08.nk'
h2O08data = np.loadtxt(h2ofile08,max_rows=486,skiprows=2)

# Call the function to get the interpolated data
new_dust_prop = interpolate_dust_properties(astrosildata, h2O08data)

plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('H2O ice')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(h2O08data[:,0], h2O08data[:,1],'C1--',label='h2O ice n')
plt.plot(h2O08data[:,0], h2O08data[:,2],'C1--',label='h2O ice k')

plt.plot(new_dust_prop[:,0], new_dust_prop[:,1],'C2--',label='h2O ice new n')
plt.plot(new_dust_prop[:,0], new_dust_prop[:,2],'C2--',label='h2O ice new k')

plt.legend()

header = "nk\nSh2O_ice 273 0.917 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/H2O_ice.nk'
np.savetxt(output_file, new_dust_prop, header=header, comments='', fmt='%.6E')

###############################################################################
# Dity Ice from https://ui.adsabs.harvard.edu/abs/1993A%26A...279..577P/abstract

dirtyfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/dirty-ice.nk'
dirtydata = np.loadtxt(dirtyfile,max_rows=90,skiprows=2)

# Call the function to get the interpolated data
new_dust_prop = interpolate_dust_properties(astrosildata, dirtydata)

plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('dirty H2O ice')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(dirtydata[:,0], dirtydata[:,1],'C1--',label='dirty h2O ice n')
plt.plot(dirtydata[:,0], dirtydata[:,2],'C1--',label='dirty h2O ice k')

plt.plot(new_dust_prop[:,0], new_dust_prop[:,1],'C2--',label='dirty h2O ice new n')
plt.plot(new_dust_prop[:,0], new_dust_prop[:,2],'C2--',label='dirty h2O ice new k')

plt.legend()


header = "nk\nSh2O_dirtyice 175 1.6 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/H2O_dirtyice.nk'
np.savetxt(output_file, new_dust_prop, header=header, comments='', fmt='%.6E')

###############################################################################
# Optical constants for forsterite (Mg2SiO4) - crystalline
# Sogawa et al. 2006, Astron. Astrophys. 451, 357

forstp1file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp1-kyoto-uvx.nk'
forstp2file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp2-kyoto-uvx.nk'
forstp3file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp3-kyoto-uvx.nk'

forstp1data = np.loadtxt(forstp1file,max_rows=6769,skiprows=2)
forstp2data = np.loadtxt(forstp2file,max_rows=6783,skiprows=2)
forstp3data = np.loadtxt(forstp3file,max_rows=6814,skiprows=2)


# Call the function to get the interpolated data
fost1 = interpolate_dust_properties(astrosildata, forstp1data)
fost2 = interpolate_dust_properties(astrosildata, forstp2data)
fost3 = interpolate_dust_properties(astrosildata, forstp3data)

# combine polarizations
wa, wb, wc = 0.6,0.2,0.2
n_eff, k_eff = get_effective_optical_constants(fost1[:,1], fost1[:,2], wa,
                                fost2[:,1], fost2[:,2], wb,
                                fost3[:,1], fost3[:,2], wc)

forstdata = np.vstack((fost1[:,0], n_eff, k_eff)).T

# plot results
plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('Fosterite')
plt.xlabel('wavelenth (microns)')

# plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n',alpha=0.5)
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k',alpha=0.5)

# plt.plot(forstp1data[:,0], forstp1data[:,1],'C1--',label='Fost. crystal n',alpha=0.5)
plt.plot(forstp1data[:,0], forstp1data[:,2],'C1--',label='Fost. crystal k',alpha=0.5)

# plt.plot(forstdata[:,0], forstdata[:,1],'C2--',label='new n')
plt.plot(forstdata[:,0], forstdata[:,2],'C2--',label='new k')

plt.legend()


header = 'nk\nSforst-cryst-kyoto-uvx 1000 3.27 0.588 20.10'

# combined file
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite_{wa}-{wb}-{wc}.nk'
np.savetxt(output_file, forstdata, header=header, comments='', fmt='%.6E')

output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite_p1.nk'
np.savetxt(output_file, fost1, header=header, comments='', fmt='%.6E')
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite_p2.nk'
np.savetxt(output_file, fost2, header=header, comments='', fmt='%.6E')
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite_p3.nk'
np.savetxt(output_file, fost3, header=header, comments='', fmt='%.6E')

###############################################################################
# Optical constants for enstatite (MgSiO3) - crystalline
# Jager et al. 1998, Astron. Astrophys. 339, 904

# enstap1file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-par-jena-uvx.nk'
# enstap2file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-perp1-jena-uvx.nk'
# enstap3file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-perp2-jena-uvx.nk'
# enstap1data = np.loadtxt(enstap1file,max_rows=1651,skiprows=2)
# enstap2data = np.loadtxt(enstap2file,max_rows=1652,skiprows=2)
# enstap3data = np.loadtxt(enstap3file,max_rows=1652,skiprows=2)
# header = 'nk\nSenst-cryst-jena-uvx 1000 3.2 0.588 20.08'
# footer=f"""
# ---------------------------------------------------------------
# weights used for each axis: {wa}, {wb}, {wc}
#  Optical constants for enstatite (MgSiO3) - crystalline.
#  Data from Jager et al. 1998, Astron. Astrophys. 339, 904 and
#    the UV-extension from Draine, 2003, ApJ 598, 1026
#  Range 0.0001-10000 microns, 1652 entries.             [NW, Apr'07]
# ---------------------------------------------------------------"""

enstap1file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/Enstatite100K_Epx.nk'
enstap2file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/Enstatite100K_Epy.nk'
enstap3file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/Enstatite100K_Epz.nk'
enstap1data = np.loadtxt(enstap1file,max_rows=5331,skiprows=2)
enstap2data = np.loadtxt(enstap2file,max_rows=5331,skiprows=2)
enstap3data = np.loadtxt(enstap3file,max_rows=5331,skiprows=2)
header = 'nk\nSEnstatite100k  1300. 3.2  0.43 25.1'
footer=f"""
---------------------------------------------------------------
weights used for each axis: {wa}, {wb}, {wc}
 Optical constants for Orthoenstatite ( Mg0.92Fe0.09SiO3) - crystalline.
 Data for cryst. silicate : https://www.astro.uni-jena.de/Laboratory/OCDB/crsilicates.html
 Ortho-enstatite (Burma, T-dependent data): (Zeidler et al. 2015)
 the UV-extension from Draine, 2003, ApJ 598, 1026
 Range 0.0001-10000 microns, 1652 entries.             [NW, Apr'07]
---------------------------------------------------------------"""



# Call the function to get the interpolated data
ensta1 = interpolate_dust_properties(astrosildata, enstap1data)
ensta2 = interpolate_dust_properties(astrosildata, enstap2data)
ensta3 = interpolate_dust_properties(astrosildata, enstap3data)

# combine polarizations
wa, wb, wc = 0.3,0.7,0.
n_eff, k_eff = get_effective_optical_constants(ensta1[:,1], ensta1[:,2], wa,
                                ensta2[:,1], ensta2[:,2], wb,
                                ensta3[:,1], ensta3[:,2], wc)

enstadata = np.vstack((ensta1[:,0], n_eff, k_eff)).T

# plot results
plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('Enstatite')
plt.xlabel('wavelenth (microns)')

# plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n',alpha=0.5)
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k',alpha=0.5)

# plt.plot(enstap1data[:,0], enstap1data[:,1],'C1--',label='ensta. crystal n',alpha=0.5)
plt.plot(enstap1data[:,0], enstap1data[:,1],'C1--',label='ensta. crystal k',alpha=0.5)

# plt.plot(enstadata[:,0], enstadata[:,1],'C2--',label='new n')
plt.plot(enstadata[:,0], enstadata[:,1],'C2--',label='new k')
plt.legend()

plt.figure()
plt.xscale('log')
plt.yscale('log')
# Calculate Qabs and other things
a_radius = 0.001

wave_micron, n_values, k_values = ensta1[:,0], ensta1[:,1], ensta1[:,2]
mie_cabs1, mie_csca1, cde_cabs1, cde_csca1, mie_qabs1, cde_qabs1, V, k_wavenumber_ref = process_nk(wave_micron, n_values, k_values, a_radius)
plt.plot(ensta1[:,0], mie_cabs1,'C1',label='Qabs p1 k',alpha=0.5)
plt.plot(ensta1[:,0], cde_cabs1,'C1--',label='Qabs CDE p1 k',alpha=0.5)

wave_micron, n_values, k_values = ensta2[:,0], ensta2[:,1], ensta2[:,2]
mie_cabs2, mie_csca2, cde_cabs2, cde_csca2, mie_qabs2, cde_qabs2, V, k_wavenumber_ref = process_nk(wave_micron, n_values, k_values, a_radius)
plt.plot(ensta2[:,0], mie_cabs2,'C2',label='Qabs p2 k',alpha=0.5)
plt.plot(ensta2[:,0], cde_cabs2,'C2--',label='Qabs CDE p2 k',alpha=0.5)

wave_micron, n_values, k_values = ensta3[:,0], ensta3[:,1], ensta3[:,2]
mie_cabs3, mie_csca3, cde_cabs3, cde_csca3, mie_qabs3, cde_qabs3, V, k_wavenumber_ref = process_nk(wave_micron, n_values, k_values, a_radius)
plt.plot(ensta3[:,0], mie_cabs3,'C3',label='Qabs p3 k',alpha=0.5)
plt.plot(ensta3[:,0], cde_cabs3,'C3--',label='Qabs CDE p3 k',alpha=0.5)

plt.vlines([40.3,43.2],0,5)

plt.legend()


# combined file
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/enstatite_{wa}-{wb}-{wc}.nk'
np.savetxt(output_file, enstadata, header=header, footer=footer, comments='', fmt='%.6E')

output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/enstatite_p1.nk'
np.savetxt(output_file, ensta1, header=header, footer=footer, comments='', fmt='%.6E')
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/enstatite_p2.nk'
np.savetxt(output_file, ensta2, header=header, footer=footer, comments='', fmt='%.6E')
output_file = f'/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/enstatite_p3.nk'
np.savetxt(output_file, ensta3, header=header, footer=footer, comments='', fmt='%.6E')


