#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 15:22:44 2023

tools to deal with the input and output of the photoionization code Mocassin

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import subprocess
import shutil

from scipy.ndimage import rotate
from scipy.spatial.distance import cdist
from numpy.lib.recfunctions import unstructured_to_structured

###############################################################################
# function to read dust species file

def read_dust_species_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    num_species = int(lines[0].strip())
    species_data = []

    for line in lines[1:num_species+1]:
        parts = line.strip().split()
        name = parts[0].strip("'")
        abundance = float(parts[1])
        rmin = float(parts[2])
        rmax = float(parts[3])
        species_data.append({
            'name': name,
            'abundance': abundance,
            'rmin': rmin,
            'rmax': rmax
        })

    return num_species, species_data

###############################################################################
# function to write dust species file

def write_dust_species_file(filepath, num_species, species_data):
    with open(filepath, 'w') as f:
        f.write(f"{num_species}\n")
        for sp in species_data:
            line = f"'{sp['name']}' {sp['abundance']}  {sp['rmin']} {sp['rmax']}\n"
            f.write(line)

###############################################################################
# function to project grid values on an image
from numba import njit, prange

@njit(parallel=True)
def splat_numba(proj_x, proj_y, intensity, volumes, img_grid_x, img_grid_y, scale_intensity):
    image = np.zeros_like(img_grid_x)
    weight_map = np.zeros_like(img_grid_x)
    cell_psf = volumes**(1/3)
    Amp = intensity * volumes if scale_intensity else intensity

    for i in prange(len(proj_x)):
        for j in range(img_grid_x.shape[0]):
            for k in range(img_grid_x.shape[1]):
                dx = img_grid_x[j, k] - proj_x[i]
                dy = img_grid_y[j, k] - proj_y[i]
                exponent = -0.5 * (dx*dx + dy*dy) / (cell_psf[i]**2)
                weight = np.exp(exponent)
                image[j, k] += Amp[i] * weight
                weight_map[j, k] += weight

    return image, weight_map

def render_view_fast(proj_x, proj_y, intensity, volumes, image_size=(50, 50), scale_intensity=False):
    norm_cst = np.max(np.abs([proj_x, proj_y]))
    proj_x, proj_y, volumes = proj_x / norm_cst, proj_y / norm_cst, volumes / norm_cst**3
    img_x = np.linspace(proj_x.min(), proj_x.max(), image_size[0])
    img_y = np.linspace(proj_y.min(), proj_y.max(), image_size[1])
    img_grid_x, img_grid_y = np.meshgrid(img_x, img_y)

    image, weight_map = splat_numba(proj_x, proj_y, intensity, volumes, img_grid_x, img_grid_y, scale_intensity)
    Amp = intensity * volumes if scale_intensity else intensity
    normalized_image = np.where(weight_map > 0, image / weight_map, 0)
    normalized_image = normalized_image / normalized_image.sum() * Amp.sum()
    return normalized_image

###############################################################################
# function to project grid values on an image

def render_view(proj_x, proj_y, intensity, volumes, image_size=50, scale_intensity=False):
    """
    Render a 2D depth-aware image from projected points.

    Parameters:
        proj_x, proj_y: 2D projected coordinates (arrays).
        intensity: Intensity values associated with each point.
        volumes: Volume values associated with each point.
        image_size: Size of the output image grid.
    """
    #normalize
    norm_cst = np.max(np.abs([proj_x, proj_y]))
    proj_x, proj_y, volumes = proj_x/norm_cst, proj_y/norm_cst, volumes/norm_cst**3
    print(norm_cst, proj_x.max(), proj_y.max(), volumes.mean())

    # Define 2D grid for binning
    img_x = np.linspace(proj_x.min(), proj_x.max(), image_size)
    img_y = np.linspace(proj_y.min(), proj_y.max(), image_size)
    img_grid_x, img_grid_y = np.meshgrid(img_x, img_y)

    cell_psf = volumes**(1/3)

    # Initialize intensity map and weight map
    image = np.zeros_like(img_grid_x)
    weight_map = np.zeros_like(img_grid_x)

    if scale_intensity:
        Amp = intensity * volumes
    else:
        Amp = intensity

    for i in range(len(proj_x)):
        exponent = -0.5 * ((img_grid_x - proj_x[i])**2 + (img_grid_y - proj_y[i])**2) / cell_psf[i]**2
        contribution = Amp[i] * np.exp(exponent)
        weight = np.exp(exponent)  # Use the same Gaussian for weight

        image = image + contribution
        weight_map = weight_map + weight

    # Normalize the image by the sum of weights
    # Avoid division by zero
    normalized_image = np.where(weight_map > 0, image / weight_map, 0)

    # conserve energy information
    normalized_image = normalized_image/normalized_image.sum() * Amp.sum()

    return normalized_image    
    
############################################################
# Function to calculate dust grain size distribution

def get_distribution(nsizes, amin, amax, slope):
    """create grain size distribution to be used in Mocassin.
       assumes correct input of values so it can be used in 
       optimization proceedure

    Args:
    Parameters:
    - nsizes (int or float): The number of size intervals to generate (will be rounded).
    - amin (float): The minimum size value (must be > 0).
    - amax (float): The maximum size value (must be >= amin).
    - slope (float): The exponent to apply to the size in the power-law distribution.

    Returns:
    - str: A formatted string showing each size and its corresponding power-law value.
           Includes a header with the number of sizes and a footer with the slope used.
           Each line in the main body contains:
             index, size (log-spaced), and size^slope in scientific notation.
    
    Notes:
    - WARNING: DOES NOT PERFORM ERROR CHECKING
       """
    try:
        nsizes = int(round(float(nsizes)))
        amin = float(amin)
        amax = float(amax)
        slope = float(slope)
    except ValueError:
        return None  # Invalid input
    

    if nsizes > 1:
        astep = (np.log10(amax) - np.log10(amin)) / (nsizes - 1)
    else:
        astep = np.log10(amax) - np.log10(amin)

    output = f"{nsizes} sizes\n\n"
    
    for i in range(nsizes):
        size = 10 ** (np.log10(amin) + i * astep)
        dist = size ** slope
        output += f"{i}\t{size:.6f}\t{dist:.6e}\n"

    output += f"\n slope (p): {slope}"
    return output


############################################################
# coordinate conversions

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

############################################################
# function to read dust temperatures from model

def read_dust_grid(filename, n_size, n_species):
    """Reads the dust grid data from a file.

    Args:
        filename: The name of the file to read.
        n_size: The number of grain size bins.
        n_species: The number of grain species.

    Returns:
        A list for:
            - 'Ndust': The dust number density for each cell.
            - 'temps': A list of grain temperatures.
    """

    temps = []
    Ndust = []
    with open(filename, 'r') as f:
        while True:
            # Read Ndust
            line = f.readline().strip()
            if not line:
                break
            Ndust.append(float(line))

            # Read grain temperatures
            grain_temps = []
            for _ in range(n_size + 1):
                line = f.readline().strip()
                grain_temps.append([float(x) for x in line.split()])

            grain_temps = np.array(grain_temps)
            temps.append(grain_temps)

    return Ndust, temps
############################################################
# make a cube from 1/8 model run with symmetricXYZ

def make_cube(cube):
    """
    Constructs a full 3D cube from a single octant, assuming the planes 
    at index 0 of the input cube are the shared central surfaces.
    """
    
    # 1. Expand along the X-axis (axis 0)
    # Flip the cube, drop the shared origin plane ([:-1, :, :]), and concatenate
    half_x = np.flip(cube, axis=0)[:-1, :, :]
    full_x = np.concatenate((half_x, cube), axis=0)
    
    # 2. Expand along the Y-axis (axis 1) using the newly expanded X-array
    # Drop the shared Y origin plane ([:, :-1, :])
    half_y = np.flip(full_x, axis=1)[:, :-1, :]
    full_xy = np.concatenate((half_y, full_x), axis=1)
    
    # 3. Expand along the Z-axis (axis 2) using the XY-array
    # Drop the shared Z origin plane ([:, :, :-1])
    half_z = np.flip(full_xy, axis=2)[:, :, :-1]
    full_xyz = np.concatenate((half_z, full_xy), axis=2)
    
    return full_xyz

############################################################
# rotate a cube aound z axis, then y axis and then x axis

def cube_rot(data,xangle,yangle,zangle,order=1):
    
    cube1 = rotate(data,zangle, axes=(0,1),prefilter=False, reshape=True)
    cube2 = rotate(cube1,yangle, axes=(1,2),prefilter=False, reshape=True)
    cube3 = rotate(cube2,xangle, axes=(0,2),prefilter=False, reshape=True)
    
    return cube3

############################################################
# Read abundance file (abundance, element)

def read_abundance_file(model_path,filename):
    
    abundances=np.genfromtxt(model_path+"input/"+filename,dtype=[('abundance', '<f8'), 
                                                                 ('element', '<U2')], 
                             encoding=None, usecols=(0,2))
    
    return abundances

############################################################
# write abundance file from recarray (abundance, element)

def write_abundance_file(model_path,filename,abund_data):
    
    np.savetxt(model_path+"input/"+filename, abund_data, fmt='%8.2e ! %s')

############################################################
# Read lineflux.out file with line fluxes
# returns list of lines split by spaces

def read_line_file(path_model):
    linefluxout = open(path_model+'output/lineFlux.out', 'r')
    lines = linefluxout.readlines()
    line_fluxes = []
    for line in lines:
        if (line.split() != []):
            line_fluxes.append(line.split())
    linefluxout.close()   
    
    return line_fluxes

############################################################
# read HI recombination lines (lambda,flux/hb)

def read_HI_lines(path_model):
    line_fluxes = read_line_file(path_model)
    inds = [i for i, val in enumerate(line_fluxes) if val[0] =='1']
    ind_start = [i for i, val in enumerate(line_fluxes) if (val[0] =='1') and 
                 (val[2] =='3') and (val[3] =='2')][0]
    ind_end = [i for i, val in enumerate(line_fluxes) if (val[0] =='1') and 
               (val[2] =='30') and (val[3] =='8')][0]
    HI_lines = np.array(line_fluxes[ind_start:ind_end+1]).astype(float)
    
    return HI_lines[:,[0,1,6,4]]

############################################################
# read HeI recombination lines (lambda,flux/hb)

def read_HeI_lines(path_model):
    line_fluxes = read_line_file(path_model)
    ind_start = [i for i, val in enumerate(line_fluxes) if (val[0] =='2') and 
                 (val[1] =='2') and (val[2] =='1')][0]
    ind_end = [i for i, val in enumerate(line_fluxes) if (val[0] =='2') and 
               (val[1] =='2') and (val[2] =='34')][0]
    HeI_lines = np.array(line_fluxes[ind_start:ind_end+1]).astype(float)
    
    return HeI_lines[:,[0,1,3,4]]

############################################################
# read HeII recombination lines (lambda,flux/hb)

def read_HeII_lines(path_model):
    line_fluxes = read_line_file(path_model)
    ind_start = [i for i, val in enumerate(line_fluxes) if (val[0:4] == ['2','3','3','2'])][0]
    ind_end = [i for i, val in enumerate(line_fluxes) if (val[0:4] == ['2','3','30','16'])][0]
    HeII_lines = np.array(line_fluxes[ind_start:ind_end+1]).astype(float)
    
    return HeII_lines[:,[0,1,6,4]]

############################################################
# read Heavy ion lines (element,lambda,flux/hb)

def read_heavyion_lines(path_model):
    line_fluxes = read_line_file(path_model)
    ind_heavy = [i for i, val in enumerate(line_fluxes) if 'Heavy' in val][0] + 1
    ind_heavy_end = [i for i, val in enumerate(line_fluxes) if 'Component:' in val or 'Recombination' in val][0] - 1
    Heavy_ion_lines = np.array(line_fluxes[ind_heavy:ind_heavy_end]).astype(float)
    
    return Heavy_ion_lines[:,[0,1,4,5]]

############################################################
# Read observations of line 
# format should be (lambda, I/Hbeta) with Hb=1

def read_obs_lines(path):
    obs_lines = np.genfromtxt(path,names=True, comments='#',dtype=None)
    return obs_lines

############################################################
# Read model lines and return a single array

def read_model_lines(path_model):
    
    hi_lines = read_HI_lines(path_model)
    hei_lines = read_HeI_lines(path_model)
    heii_lines = read_HeII_lines(path_model)
    heavy_lines = read_heavyion_lines(path_model)
    
    model_lines = np.concatenate([hi_lines,hei_lines,heii_lines,heavy_lines])
    
    dt = np.dtype([('A', '<i8'), ('lowerLevel', '<i8'), ('lambda', '<f8'), ('intensity', '<f8')])
    model_lines = unstructured_to_structured(model_lines, dtype=dt)
    
    return model_lines
    
############################################################
# Cross match model to observational line intensities
# !!!!WARNING!!!! -  DOES NOT HANDLE BLENDS (YET)

def Xmatch_lines(obs_data,mod_data,tol=3.0, blends=[]):
    
    line_data = []
        
    for i,line in enumerate(obs_data['lambda']):
        
        ind = np.argmin(np.abs(line-mod_data['lambda']) + np.abs(obs_data['A'][i]-mod_data['A']))
        dl = np.abs(line-mod_data['lambda'][ind]) + np.abs(obs_data['A'][i]-mod_data['A'][ind])
        
        if (dl < tol):
                        
            if line in blends:
                #print(dl,mod_data['lambda'][np.abs(line-mod_data['lambda']) < tol])
                int_sum = np.sum(mod_data['intensity'][np.abs(line-mod_data['lambda']) < tol])
                lambda_mean = np.mean(mod_data['lambda'][np.abs(line-mod_data['lambda']) < tol])
                                                       
                line_data.append([obs_data['lambda'][i], lambda_mean, 
                                  obs_data['intensity'][i], obs_data['error'][i], 
                                  obs_data['deltaO'][i], int_sum])
            else:
                line_data.append([obs_data['lambda'][i], mod_data['lambda'][ind], 
                                  obs_data['intensity'][i], obs_data['error'][i],
                                  obs_data['deltaO'][i], mod_data['intensity'][ind]])
                
                

            # with np.printoptions(precision=6, suppress=True):
            #     print(line, mod_data[ind])
                
    line_data = np.array(line_data)
    dt = np.dtype([('lambda', '<f8'), ('rest', '<f8'), ('I_obs', '<f8'), 
                   ('I_obs_er', '<f8'), ('I_obs_dO', '<f8'), ('I_mod', '<f8')])
    line_data = unstructured_to_structured(line_data, dtype=dt)
    
    return line_data


############################################################
# Run mocassin

def run_mocassin(path_model,ncpu=1,verbose=False):
    """
    Executes the MOCASSIN photoionization code using MPI parallelization.

    Parameters:
    - path_model (str): The path to the model directory containing the 'input/' subdirectory 
                        and where the output file 'saida.dat' will be written.
    - ncpu (int, optional): Number of CPUs to use in the MPI execution. Default is 1.
    - verbose (bool, optional): If True, enables verbose output (currently not used).

    Raises:
    - FileNotFoundError: If the required 'input/' subdirectory does not exist inside path_model.
    - subprocess.CalledProcessError: If the MOCASSIN execution fails.

    Behavior:
    - Verifies the existence of the 'input/' directory in the specified model path.
    - Checks if the `mpirun` command is available in the system's PATH.
    - Executes MOCASSIN using `mpirun`, passing the number of CPUs and expecting a hostfile named "machines".
    - Captures standard output from the simulation and writes it to 'saida.dat' in the model directory.

    Notes:
    - The path to the mocassin binary is hardcoded: '/home/hmonteiro/mocassin_2023/bin/mocassin'.
      Adjust this if your binary is located elsewhere.
    """
    
    # check to see if input directory exists
    if not os.path.exists(path_model+'input/') or not os.path.isdir(path_model+'input/'):
        raise FileNotFoundError(f"The directory (input) with parameter definitions does not exist.")
        
    # Check if mpirun exists in the system's PATH
    if shutil.which('mpirun') is not None:
        mpirun_path = shutil.which('mpirun')
        # print(f"The program mpirun exists in the system.")
    else:
        print(f"The program mpirun does not exist in the system.")    

    # Run the code
    ncpu = str(ncpu)
    command = ["mpirun", "-np", str(ncpu), "--hostfile", "machines", "/home/hmonteiro/mocassin_2023/bin/mocassin"]

    try:
        with open(path_model+'saida.dat', 'w') as fd:
            subprocess.run(command, check=True, cwd=path_model, stdout=fd)
        # print("mocassin execution successful.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing mocassin: {e}")
    

###############################################################################
# likelihood for emission line fit  with errors. The variable obs_data is read 
# with 'read_obs_lines'

def log_likelihood(theta, obs_data, weight, mod_dir, use_kappa=True):
    
    
    He_abun, C_abun, N_abun, O_abun, S_abun = theta

    ###########################################################################
    # read abundnce file
    abun_file = 'abun-ngc3132.in'
    mod_abund = read_abundance_file(mod_dir,abun_file)
    
    #change abundance values    
    ind = np.where(mod_abund['element'] == 'He')
    mod_abund['abundance'][ind] = He_abun
    
    ind = np.where(mod_abund['element'] == 'C')
    mod_abund['abundance'][ind] = C_abun*1.0e-4
    
    ind = np.where(mod_abund['element'] == 'N')
    mod_abund['abundance'][ind] = N_abun*1.0e-4
    
    ind = np.where(mod_abund['element'] == 'O')
    mod_abund['abundance'][ind] = O_abun*1.0e-4
    
    ind = np.where(mod_abund['element'] == 'S')
    mod_abund['abundance'][ind] = S_abun*1.0e-5
    
    # write a modified file
    write_abundance_file(mod_dir, abun_file, mod_abund)
    
    ###########################################################################    
    # run mocassin
    run_mocassin(mod_dir,ncpu=120)
    
    ###########################################################################    
    # read model results
    mod_lines = read_model_lines(mod_dir)
    line_data = Xmatch_lines(obs_data,mod_lines,tol=3.0, blends=[1907.,3727.,7322.,7332.,5200])

    ###########################################################################    
    # use diagnostic ratios
    
    if (6718. in line_data['lambda']):
        
        ind_6718 = np.where(line_data['lambda'] == 6718.)[0][0] 
        I_6718 = line_data[ind_6718]['I_obs']
        
        ind_6732 = np.where(line_data['lambda'] == 6732.)[0][0] 
        I_6732 = line_data[ind_6732]['I_obs']
        
        s_diag_obs = I_6732/I_6718
        s_diag_obs_er = s_diag_obs*np.sqrt((line_data[ind_6718]['I_obs_er']/line_data[ind_6718]['I_obs'])**2 + 
                                (line_data[ind_6732]['I_obs_er']/line_data[ind_6732]['I_obs'])**2)
        
        I_6718_mod = line_data[ind_6718]['I_mod']
        I_6732_mod = line_data[ind_6732]['I_mod']
        s_diag_mod = I_6732_mod/I_6718_mod
                
    if (4364. in line_data['lambda']):  # (4959+5007)/4363
        ind_4364 = np.where(line_data['lambda'] == 4364.)[0][0] 
        I_4364 = line_data[ind_4364]['I_obs']
        
        ind_5008 = np.where(line_data['lambda'] == 5008.)[0][0] 
        I_5008 = line_data[ind_5008]['I_obs']
        
        o3_diag_obs = (4/3*I_5008)/I_4364 
        o3_diag_obs_er = o3_diag_obs*np.sqrt((line_data[ind_4364]['I_obs_er']/line_data[ind_4364]['I_obs'])**2 +
                                             (line_data[ind_5008]['I_obs_er']/line_data[ind_5008]['I_obs'])**2)
        
        I_4364_mod = line_data[ind_4364]['I_mod']
        I_5008_mod = line_data[ind_5008]['I_mod']
        o3_diag_mod = (4/3*I_5008_mod)/I_4364_mod 
        
    if (5756. in line_data['lambda']):  # (6584+6546)/5754
        ind_5756 = np.where(line_data['lambda'] == 5756.)[0][0] 
        I_5756 = line_data[ind_5756]['I_obs']
        
        ind_6585 = np.where(line_data['lambda'] == 6585.)[0][0] 
        I_6585 = line_data[ind_6585]['I_obs']
        
        n2_diag_obs = (4/3*I_6585)/I_5756 
        n2_diag_obs_er = n2_diag_obs*np.sqrt((line_data[ind_5756]['I_obs_er']/line_data[ind_5756]['I_obs'])**2 + 
                                 (line_data[ind_6585]['I_obs_er']/line_data[ind_6585]['I_obs'])**2)
        
        I_5756_mod = line_data[ind_5756]['I_mod']
        I_6585_mod = line_data[ind_6585]['I_mod']
        n2_diag_mod = (4/3*I_6585_mod)/I_5756_mod 
        
    if (6312. in line_data['lambda']):  # [SIII](9069+9531)/6312
        ind_6312 = np.where(line_data['lambda'] == 6312.)[0][0] 
        I_6312 = line_data[ind_6312]['I_obs']
        
        ind_9072 = np.where(line_data['lambda'] == 9072.)[0][0] 
        I_9072 = line_data[ind_9072]['I_obs']
        
        s3_diag_obs = (3.48*I_9072)/I_6312 
        s3_diag_obs_er = s3_diag_obs*np.sqrt((line_data[ind_6312]['I_obs_er']/line_data[ind_6312]['I_obs'])**2 + 
                                 (line_data[ind_9072]['I_obs_er']/line_data[ind_9072]['I_obs'])**2)
        
        I_6312_mod = line_data[ind_6312]['I_mod']
        I_9072_mod = line_data[ind_9072]['I_mod']
        s3_diag_mod = (3.48*I_9072_mod)/I_6312_mod 
        
        
    Diags = [(s_diag_mod-s_diag_obs)**2/s_diag_obs_er**2,
             (o3_diag_mod-o3_diag_obs)**2/o3_diag_obs_er**2,
             (n2_diag_mod-n2_diag_obs)**2/n2_diag_obs_er**2,
             (s3_diag_mod-s3_diag_obs)**2/s3_diag_obs_er**2]
    
    #print(Diags)
    #print(s_diag_obs, o3_diag_obs, n2_diag_obs)
    #print(s_diag_obs_er,o3_diag_obs_er,n2_diag_obs_er)

    if use_kappa:
        
        # use metric recomended by https://arxiv.org/abs/2312.01873
        ###########################################################################    
        # Calculate likelihood using kappa metric
        
        chi = (line_data['I_obs']-line_data['I_mod'])**2/line_data['I_obs_er']**2    
        tau = np.log10(1+line_data['I_obs_dO']/line_data['I_obs'])
        kappa = (np.log10(line_data['I_mod']) - np.log10(line_data['I_obs']))/tau  
        
    else:
        ###########################################################################    
        # Calculate likelihood using chi2     
            
        S = (line_data['I_mod'] - line_data['I_obs'])**2 / line_data['I_obs_er']**2
        rel_er = (line_data['I_mod'] - line_data['I_obs'])**2/line_data['I_obs']**2    
        logL = -0.5 * np.sum(S)   
        logL_diags = -0.5*np.sum(Diags)


        
    return - (logL + logL_diags)




###############################################################################
# Function to print out table of results and goodness of fit

def print_tab_results(obs_file,mod_dir):
    
    ###########################################################################
    # Read observed line intensities
    obs_data = read_obs_lines(obs_file)   
        
    ###########################################################################    
    # read model results
    mod_lines = read_model_lines(mod_dir)
    line_data = Xmatch_lines(obs_data,mod_lines,tol=3.0, blends=[1907.,3727.,7322.,7332.,5200])

    ###########################################################################    
    # use diagnostic ratios
    
    if (6718. in line_data['lambda']):
        
        ind_6718 = np.where(line_data['lambda'] == 6718.)[0][0] 
        I_6718 = line_data[ind_6718]['I_obs']
        
        ind_6732 = np.where(line_data['lambda'] == 6732.)[0][0] 
        I_6732 = line_data[ind_6732]['I_obs']
        
        s_diag_obs = I_6732/I_6718
        s_diag_obs_er = s_diag_obs*np.sqrt((line_data[ind_6718]['I_obs_er']/line_data[ind_6718]['I_obs'])**2 + 
                                (line_data[ind_6732]['I_obs_er']/line_data[ind_6732]['I_obs'])**2)
        
        I_6718_mod = line_data[ind_6718]['I_mod']
        I_6732_mod = line_data[ind_6732]['I_mod']
        s_diag_mod = I_6732_mod/I_6718_mod
                
    if (4364. in line_data['lambda']):  # (4959+5007)/4363
        ind_4364 = np.where(line_data['lambda'] == 4364.)[0][0] 
        I_4364 = line_data[ind_4364]['I_obs']
        
        ind_5008 = np.where(line_data['lambda'] == 5008.)[0][0] 
        I_5008 = line_data[ind_5008]['I_obs']
        
        o3_diag_obs = (4/3*I_5008)/I_4364 
        o3_diag_obs_er = o3_diag_obs*np.sqrt((line_data[ind_4364]['I_obs_er']/line_data[ind_4364]['I_obs'])**2 +
                                             (line_data[ind_5008]['I_obs_er']/line_data[ind_5008]['I_obs'])**2)
        
        I_4364_mod = line_data[ind_4364]['I_mod']
        I_5008_mod = line_data[ind_5008]['I_mod']
        o3_diag_mod = (4/3*I_5008_mod)/I_4364_mod 
        
    if (5756. in line_data['lambda']):  # (6584+6546)/5754
        ind_5756 = np.where(line_data['lambda'] == 5756.)[0][0] 
        I_5756 = line_data[ind_5756]['I_obs']
        
        ind_6585 = np.where(line_data['lambda'] == 6585.)[0][0] 
        I_6585 = line_data[ind_6585]['I_obs']
        
        n2_diag_obs = (4/3*I_6585)/I_5756 
        n2_diag_obs_er = n2_diag_obs*np.sqrt((line_data[ind_5756]['I_obs_er']/line_data[ind_5756]['I_obs'])**2 + 
                                 (line_data[ind_6585]['I_obs_er']/line_data[ind_6585]['I_obs'])**2)
        
        I_5756_mod = line_data[ind_5756]['I_mod']
        I_6585_mod = line_data[ind_6585]['I_mod']
        n2_diag_mod = (4/3*I_6585_mod)/I_5756_mod 
        
    if (6312. in line_data['lambda']):  # [SIII](9069+9531)/6312
        ind_6312 = np.where(line_data['lambda'] == 6312.)[0][0] 
        I_6312 = line_data[ind_6312]['I_obs']
        
        ind_9072 = np.where(line_data['lambda'] == 9072.)[0][0] 
        I_9072 = line_data[ind_9072]['I_obs']
        
        s3_diag_obs = (3.48*I_9072)/I_6312 
        s3_diag_obs_er = s3_diag_obs*np.sqrt((line_data[ind_6312]['I_obs_er']/line_data[ind_6312]['I_obs'])**2 + 
                                 (line_data[ind_9072]['I_obs_er']/line_data[ind_9072]['I_obs'])**2)
        
        I_6312_mod = line_data[ind_6312]['I_mod']
        I_9072_mod = line_data[ind_9072]['I_mod']
        s3_diag_mod = (3.48*I_9072_mod)/I_6312_mod 
        
        
    chi_diags = [(s_diag_mod-s_diag_obs)**2/s_diag_obs_er**2,
             (o3_diag_mod-o3_diag_obs)**2/o3_diag_obs_er**2,
             (n2_diag_mod-n2_diag_obs)**2/n2_diag_obs_er**2,
             (s3_diag_mod-s3_diag_obs)**2/s3_diag_obs_er**2]
    
    ###########################################################################    
    # Calculate likelihood using chi2     
        
    S = (line_data['I_mod'] - line_data['I_obs'])**2 / line_data['I_obs_er']**2
    rel_er = (line_data['I_mod'] - line_data['I_obs'])**2/line_data['I_obs']**2    
    logL = -0.5 * np.sum(S)
    logL_diags = -0.5*np.sum(chi_diags)
        
    ###########################################################################    

    chi = (line_data['I_obs']-line_data['I_mod'])**2/line_data['I_obs_er']**2 
    dOoO = line_data['I_obs_dO']/line_data['I_obs']
    tau = np.log10(1+dOoO)
    kappa = (np.log10(line_data['I_mod']) - np.log10(line_data['I_obs']))/tau    
    rel_err = (line_data['I_obs']-line_data['I_mod'])/line_data['I_obs']

    ###########################################################################
    ## Including Chi    
    # print('')
    # print('  Line       OBS       Mod    rel. error    chi    tau')
    # print('------------------------------------------------')
    
    # for i in range(line_data.size):
    #     print('%8.1f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f'%(line_data['lambda'][i], 
    #                                                line_data['I_obs'][i],line_data['I_mod'][i], rel_err[i],chi[i], kappa[i]))

    ###########################################################################    
    print('')
    print('  Line       OBS       Error     dO/O      Model     Kappa')
    print('------------------------------------------------')
    
    for i in range(line_data.size):
        print('%8.1f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f'%(line_data['lambda'][i], 
                                                   line_data['I_obs'][i],line_data['I_obs_er'][i], dOoO[i],
                                                   line_data['I_mod'][i], kappa[i]))

    ###########################################################################    
    # for diagnostic ratios
    diags_label = ['[SII]6731/6717','[OIII](4959+5007)/4363','[NII](6584+6546)/5754','[SIII](9069+9531)/6312']
    obs_diags = np.array([s_diag_obs, o3_diag_obs, n2_diag_obs, s3_diag_obs])
    obs_diags_er = np.array([s_diag_obs_er, o3_diag_obs_er, n2_diag_obs_er, s3_diag_obs_er])
    mod_diags = np.array([s_diag_mod, o3_diag_mod, n2_diag_mod, s3_diag_mod])
    
    dOoO = 3*obs_diags_er/obs_diags
    tau = np.log10(1+dOoO)
    kappa = (np.log10(mod_diags) - np.log10(obs_diags))/tau    
    
    print('------------------------------------------------------------------')
    print('Line Diagnostics')
    for i in range(len(diags_label)):
        print('%-25s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%(diags_label[i], obs_diags[i], obs_diags_er[i], dOoO[i], mod_diags[i],kappa[i]) )
                      

    return 
















