#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:02:18 2025

@author: hmonteiro

NOTES:
    
size parameter defined in Bohren and Huffman (1983) - x = k*a = 2*pi*N*a/lambda -> a is the radius of the sphere
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

###############################################################################

def BHmie(x, refrel):
    # Ensure double precision
    #x = np.float64(x_in)
    #refrel = np.complex128(refrel_in)

    nang = 2
    dx = x
    y = x * refrel

    xstop = x + 4.0 * x**(1.0 / 3.0) + 2.0
    nstop = int(xstop)
    ymod = abs(y)
    nmx = int(max(xstop, ymod)) + 15

    dang = np.pi / 2.0 / (nang - 1)

    # 1-based indexing: index 0 is unused
    theta = np.zeros(nang + 1)
    amu = np.zeros(nang + 1)
    for j in range(1, nang + 1):
        theta[j] = (j - 1) * dang
        amu[j] = np.cos(theta[j])

    d = np.zeros(nmx + 1, dtype=complex)
    d[nmx] = 0.0 + 0.0j

    for n in range(1, nmx):
        rn =nmx - n + 1
        d[nmx - n] = (rn / y) - (1.0 / (d[nmx - n + 1] + rn / y))

    pi0 = np.zeros(nang + 1)
    pi1 = np.ones(nang + 1)
    pii = np.zeros(nang + 1)
    tau = np.zeros(nang + 1)

    nn = 2 * nang - 1
    s1 = np.zeros(nn + 1, dtype=complex)
    s2 = np.zeros(nn + 1, dtype=complex)

    psi0 = np.cos(dx)
    psi1 = np.sin(dx)
    chi0 = -np.sin(x)
    chi1 = np.cos(x)
    apsi0 = psi0
    apsi1 = psi1

    xi0 = np.complex128(apsi0 - 1j * chi0)
    xi1 = np.complex128(apsi1 - 1j * chi1)

    qsca = 0.0
    ggsca = 0.0
    n = 1
    an1 = 0.0 + 0.0j
    bn1 = 0.0 + 0.0j

    while True:
        dn = n
        rn = n
        fn = (2.0 * rn + 1.0) / (rn * (rn + 1.0))

        psi = (2.0 * dn - 1.0) * psi1 / dx - psi0
        apsi = psi
        chi = (2.0 * rn - 1.0) * chi1 / x - chi0
        xi = np.complex128(apsi - 1j * chi)

        if n > 1:
            an1 = an
            bn1 = bn

        an = ((d[n] / refrel + rn / x) * apsi - apsi1) / ((d[n] / refrel + rn / x) * xi - xi1)
        bn = ((refrel * d[n] + rn / x) * apsi - apsi1) / ((refrel * d[n] + rn / x) * xi - xi1)

        qsca += np.real((2.0 * rn + 1.0) * (abs(an)**2 + abs(bn)**2))
        ggsca += ((2.0 * rn + 1.0) / (rn * (rn + 1.0))) * (
            an.real * bn.real + an.imag * bn.imag
        )

        if n > 1:
            ggsca += ((rn - 1.0) * (rn + 1.0) / rn) * (
                an1.real * an.real + an1.imag * an.imag +
                bn1.real * bn.real + bn1.imag * bn.imag
            )

        for j in range(1, nang + 1):
            jj = 2 * nang - j
            pii[j] = pi1[j]
            tau[j] = rn * amu[j] * pii[j] - (rn + 1.0) * pi0[j]
            p = (-1.0)**(n - 1)
            t = (-1.0)**n
            s1[j] += fn * (an * pii[j] + bn * tau[j])
            s2[j] += fn * (an * tau[j] + bn * pii[j])

            if j == jj:
                break

            s1[jj] += fn * (an * pii[j] * p + bn * tau[j] * t)
            s2[jj] += fn * (an * tau[j] * t + bn * pii[j] * p)

        psi0 = psi1
        psi1 = psi
        apsi1 = psi1
        chi0 = chi1
        chi1 = chi
        xi1 = np.complex128(apsi1 - 1j * chi1)

        n += 1
        rn = n
        for j in range(1, nang + 1):
            pi1[j] = ((2.0 * rn - 1.0) / (rn - 1.0)) * amu[j] * pii[j]
            pi1[j] -= rn * pi0[j] / (rn - 1.0)
            pi0[j] = pii[j]

        if n - 1 - nstop >= 0:
            break

    ggsca = 2.0 * ggsca / qsca
    qsca = 2.0 * qsca / (x * x)
    qext = 4.0 / (x * x) * s1[1].real

    return qext, qsca, ggsca


##############################################################################
    """
    From Dr. Barlows email
    
    For Rayleigh:
    Cabs = k * V * imag(3*alpha)
    Csca = k^4 * V^2 * 3 * absolute(alpha)^2 / 2*pi
    Where k is the wavenumber ( 2pi/lambda(cm) ), V is the particle volume (cm^3)
    and alpha is given by:   alpha = (m^2 - 1) / (m^2 + 2), and m=N+iK (the optical
    constants, note: different K).
    
    For CDE:
    Cabs = k * V * imag(alpha)
    Csca = k * V * imag(alpha) / sigma
    where k and V are the same as above, and alpha is now given by:
    alpha = (2*m^2)*(ln(m^2))/(m^2 - 1) - 2, and sigma is:
    sigma = (6*pi/(V*k^3))*(imag(m^2))/(absolute(m^2 - 1))^2

    """

def calculate_Rayleigh_cabs(k, V, alpha):
    """
    Calculates Mie absorption cross-section (Cabs).
    Cabs = k * V * imag(3*alpha)
    """
    return k * V * alpha.imag * 3

def calculate_Rayleigh_csca(k, V, alpha):
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

def cde_efficiency(a,refIndex, wavelength):

    # Complex refractive index 'm'
    m = refIndex
    k_wavenumber = (2 * np.pi) / wavelength

    # Assuming spherical particle volume
    V = (4/3) * np.pi * (a**3)
    # Calculate geometrical cross-section (Ag)
    Ag = np.pi * (a**2)

    alpha = (2 * m**2 * np.log(m**2)) / (m**2 - 1) - 2
    sigma_denominator = (abs(m**2 - 1))**2
    sigma = (6 * np.pi / (V * k_wavenumber**3)) * (m**2).imag / sigma_denominator
    
    
            
    cde_cabs = k_wavenumber * V * alpha.imag 
    cde_csca = k_wavenumber * V * alpha.imag / sigma 
    
    cde_qsca = cde_csca / Ag
    cde_qabs = cde_cabs / Ag 
    cde_qext = cde_qsca + cde_qabs
    
    if ((2 * np.pi * a) / wavelength > 1):
        cde_qext = 2
        cde_qsca = cde_qext - cde_qabs
    
    return cde_qsca, cde_qabs, cde_qext

###############################################################################
# Read original Mocassin calculation

mocfile = '/home/hmonteiro/Downloads/MIe_Codes/mocassin_bhmie_ori.dat'
mocdata = np.loadtxt(mocfile)



fr1Ryd = 3.28984e15         # frequency at 1 Ryd [Hz]

#particle_radius_microns = 0.913033009 # Example: 0.1 microns
particle_radius_microns = 0.438697994 # Example: 0.1 microns
#particle_radius_microns = 1. # Example: 0.1 microns


nuMax = 15.
nuMin = 1.0e-6
nbins = 700
nuArray = np.linspace(nuMin,nuMax,nbins)
nuArray = mocdata[:,1]

# size parameter
sizeParam=2.0*3.14159265*particle_radius_microns/(2.9979250e14/(nuArray*fr1Ryd) )
sizeParam[sizeParam>100]=np.float64(100)

file_path = '/home/hmonteiro/Downloads/mocassin-rw_changes_2023/dustData/sil-dlee.nk'
file_path = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite_p1.nk'

wavelengths_micron = []
n_values = []
k_values = []

with open(file_path, 'r') as f:
    lines = f.readlines()

# Skip header lines (first two lines)
data_lines = lines[2:len(lines)-8]

for line in data_lines:
    parts = line.split()
    if len(parts) == 3:
        wavelengths_micron.append(float(parts[0]))
        n_values.append(float(parts[1]))
        k_values.append(float(parts[2]))


plt.figure()
plt.plot(wavelengths_micron,n_values)
plt.plot(wavelengths_micron,k_values)
plt.yscale('log') 
plt.xscale('log')

# Interpolate optical data in the nuArray that will be used
x_ryd = 2.9979250e14/np.array(wavelengths_micron)/fr1Ryd

n_values_new = np.interp(nuArray,x_ryd[np.argsort(x_ryd)], np.array(n_values)[np.argsort(x_ryd)])
k_values_new = np.interp(nuArray,x_ryd[np.argsort(x_ryd)], np.array(k_values)[np.argsort(x_ryd)])

##############################################################################
#%%
# Convert particle radius to cm for volume calculation
particle_radius_cm = particle_radius_microns * 1e-4 # 1 micron = 1e-4 cm
# Assuming spherical particle volume
V = (4/3) * np.pi * (particle_radius_cm**3)
# Calculate geometrical cross-section (Ag)
Ag = np.pi * (particle_radius_cm**2)


# Lists to store calculated cross-sections and efficiencies
rayl_cabs_list = []
rayl_qsca_list = []
rayl_qabs_list = []

cde_cabs_list = []
cde_qsca_list = []
cde_qabs_list = [] 

bhmie_qsca_list = [] 
bhmie_qext_list = [] 

ray_qsca_list = [] 
ray_qext_list = [] 

# get nuarray as in mocassin
lamb_hz = 2.9979250e14/np.array(wavelengths_micron)  # convert wav to freq. - from microns to hertz
lamb_ryd = lamb_hz/fr1Ryd


for i in range(nuArray.size):
    lam_micron = (2.9979250e14/(nuArray*fr1Ryd))[i]
    n_val = n_values_new[i]
    k_val = k_values_new[i]

    # Convert wavelength from microns to cm for wavenumber calculation
    lam_cm = lam_micron * 1e-4 # 1 micron = 1e-4 cm
    k_wavenumber = (2 * np.pi) / lam_cm # k is wavenumber (2pi/lambda(cm))

    # Complex refractive index 'm'
    m = complex(n_val, k_val)

    # --- Mie Theory Calculations ---
    # alpha for Mie theory
    mie_alpha = (m**2 - 1) / (m**2 + 2)
    rayl_cabs = calculate_Rayleigh_cabs(k_wavenumber, V, mie_alpha)
    rayl_csca = calculate_Rayleigh_csca(k_wavenumber, V, mie_alpha)
    rayl_cabs_list.append(rayl_cabs)
    rayl_qsca_list.append(rayl_csca / Ag)
    # Calculate Mie Absorption Efficiency (Qabs)
    rayl_qabs_list.append(rayl_cabs / Ag if Ag != 0 else np.nan)
    #print(lamb_ryd[i], lam_micron, mie_cabs, mie_csca)
    
    # BHMie
    bhmie_qext,bhmie_qsca,bhmie_gsca = BHmie(sizeParam[i],complex(n_val,k_val))
    if (sizeParam[i] == 2.6741287597254515):
        print(sizeParam[i],n_val,k_val,bhmie_qext, bhmie_qsca)
    bhmie_qext_list.append(bhmie_qext)
    bhmie_qsca_list.append(bhmie_qsca)
    
    
    # --- CDE Theory Calculations ---
    # alpha for CDE theory
    # Handle log(m^2) carefully, especially if m^2 is zero or negative real.
    # cmath.log handles principal value (n=0) by default.
    try:
        cde_alpha = (2 * m**2 * np.log(m**2)) / (m**2 - 1) - 2
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

    # cde_cabs = calculate_cde_cabs(k_wavenumber, V, cde_alpha)
    # cde_csca = calculate_cde_csca(k_wavenumber, V, cde_alpha, sigma)
    # cde_cabs_list.append(cde_cabs)
    # cde_qsca_list.append(cde_csca / Ag)
    # # Calculate CDE Absorption Efficiency (Qabs)
    # cde_qabs_list.append(cde_cabs / Ag if Ag != 0 else np.nan)
    
    cde_qsca, cde_qabs, cde_qext = cde_efficiency(particle_radius_microns,m, lam_micron)
    cde_qabs_list.append(cde_qabs)
    cde_qsca_list.append(cde_qsca)

# plt.figure()
# plt.plot((2.9979250e14/(nuArray*fr1Ryd)),mie_cabs_list,label='Cabs')
# plt.plot((2.9979250e14/(nuArray*fr1Ryd)),mie_csca_list,label='Csca')

# plt.yscale('log') 
# plt.xscale('log')
# plt.legend()


sizeParam=2.0*3.14159265*particle_radius_microns/( 2.9979250e14/(np.array(nuArray)*fr1Ryd) )
sizeParam[sizeParam>100]=100
# plt.figure()
# plt.plot(wavelengths_micron, sizeParam)
# plt.yscale('log') 
# plt.xscale('log')

# # plot efficiencies

plt.figure()
# plt.plot(2.9979250e14/mocdata[:,1]/fr1Ryd, mocdata[:,2], label='Moc Qext')
#plt.plot(2.9979250e14/mocdata[:,1]/fr1Ryd, mocdata[:,3], label='Moc Qsca')

plt.plot((2.9979250e14/(nuArray*fr1Ryd)),bhmie_qext_list,'C0-',label='Mie Qext',alpha=0.9)
#plt.plot((2.9979250e14/(nuArray*fr1Ryd)),bhmie_qsca_list,'C1-',label='Mie Qsca',alpha=0.5)

# mie_qext_list = np.array(mie_qsca_list) + np.array(mie_qabs_list)
# plt.plot((2.9979250e14/(nuArray*fr1Ryd)),ray_qext_list,'-',label='Ray. Qext')
# plt.plot((2.9979250e14/(nuArray*fr1Ryd)),ray_qsca_list,'-',label='Ray Qsca')

rayl_qext_list = np.array(rayl_qsca_list) + np.array(rayl_qabs_list)
plt.plot((2.9979250e14/(nuArray*fr1Ryd)),rayl_qext_list,'C2:',label='Rayl. Qext',alpha=0.5)
#plt.plot((2.9979250e14/(nuArray*fr1Ryd)),rayl_qsca_list,'C1:',label='Rayl. Qsca')

cde_qext_list = np.array(cde_qsca_list) + np.array(cde_qabs_list)
plt.plot((2.9979250e14/(nuArray*fr1Ryd)),cde_qext_list,'C1',label='CDE Qext')
#plt.plot((2.9979250e14/(nuArray*fr1Ryd)),cde_qsca_list,'C1--',label='CDE Qsca')
#plt.plot((2.9979250e14/(nuArray*fr1Ryd)),cde_qabs_list,'C1--',label='CDE Qabs')

# plt.plot((2.9979250e14/(nuArray*fr1Ryd)),sizeParam,'C3-',label='x',alpha=0.5)
plt.axvspan((2.9979250e14/(nuArray[sizeParam<=1][-1]*fr1Ryd)), 
            (2.9979250e14/(nuArray*fr1Ryd)).max(), color='gray', alpha=0.3, label=r'$2 \pi a / \lambda < 1$')
plt.yscale('log') 
plt.xscale('log')
plt.xlim(0.01,1e3)
plt.ylim(1e-4,1e2)
plt.xlabel('Wavelength (microns)')
plt.ylabel('Efficiency (dimensionless)')
plt.legend()


