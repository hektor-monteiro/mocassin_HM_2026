#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 17:25:52 2024

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys, os

# # add mocassin tools module
# sys.path.insert(0, 'mocassin_fit')
# from mocassin_tools import *

#plt.close('all')
############################################################
def read_dust_grid(filename, n_size, n_species):
    """Reads the dust grid data from a file.

    Args:
        filename: The name of the file to read.
        n_size: The number of grain size bins.
        n_species: The number of grain species.

    Returns:
        A list of dictionaries, where each dictionary represents a cell's data.
        Each dictionary contains the following keys:
            - 'Ndust': The dust number density.
            - 'grain_temps': A 2D NumPy array of grain temperatures.
    """

    temps = []
    Ndust = []
    with open(filename, 'r') as f:
        while True:
            # Read Ndust
            line = f.readline().strip().split()
            if not line:
                break
            Ndust.append(float(line[0]))

            # Read grain temperatures
            grain_temps = []
            for _ in range(n_size + 1):
                line = f.readline().strip()
                grain_temps.append([float(x) for x in line.split()])

            grain_temps = np.array(grain_temps)
            temps.append(grain_temps)

    return Ndust, temps

############################################################
# Reading Grid Structure
        
print ('Reading Grid Structure ...')
grid0 = open("output/grid0.out", "r")
lines = grid0.readlines()

nx,ny,nz, aux = [int(x) for x in (lines[1].split())[0:4]]
ncell = nx*ny*nz

aux,Rout = [float(x) for x in (lines[1].split())[4:]]
X = np.array([float(x) for x in lines[2:2+nx]])
Y = np.array([float(x) for x in lines[2+nx:2+nx+ny]])
Z = np.array([float(x) for x in lines[2+nx+ny:2+nx+ny+nz]])

XX, YY, ZZ = np.meshgrid(X,Y,Z, indexing='ij')
Radius = np.sqrt(XX**2 + YY**2 + ZZ**2)

aux = np.array([[int(y) for y in x.split()] for x in lines[2+nx+ny+nz:]])
converged = np.reshape(aux[:,1],(nx,ny,nz))

plt.figure()
plt.imshow(converged[:,:,int(nx/2)])
plt.title('Cell convergente')

if (nx == ny == nz or Rout == 0):
    Rout = np.max([np.abs(X).max(),np.abs(Y).max(),np.abs(Z).max()]) 

############################################################
# Reading Grid cell volumes
        
print ('Reading Grid cell volumes ...')
data = np.loadtxt("output/grid4.out")
dV = np.reshape(data[:,3],(nx,ny,nz))

############################################################
#%%

ndust, temps = read_dust_grid('output/dustGrid.out',35,2)

# nx,ny,nz = 15, 15, 21
# nx,ny,nz = 39, 39, 51

Ndust = np.reshape(np.array(ndust),(nx,ny,nz))

# Temps = [cells, radii, species]
# [:,1:,1] - silicates
# [:,1:,2] - PaH

species = 2

aux=np.mean(np.array(temps)[:,1:,species],axis=1,where=np.array(temps)[:,1:,species] > 0.)
avg_Tdust = np.reshape(aux,(nx,ny,nz))
avg_Tdust[~np.isfinite(avg_Tdust)]=0

levels=np.logspace(np.log10(avg_Tdust.max()*0.01),np.log10(avg_Tdust.max()*0.99),50)


f,ax = plt.subplots(1,3,figsize=(17,4))
p1 = ax[0].contourf(Z,X,avg_Tdust[nx//2,:,:],levels=levels, norm=LogNorm())
cb=f.colorbar(p1, ax=ax[0])
cb.set_label('T dust (k)')
ax[0].set_title('X axis cut')
ax[0].axis('off')
ax[0].set_aspect(1)
p2 = ax[1].contourf(Z,X,avg_Tdust[:,ny//2,:],levels=levels, norm=LogNorm())
cb=f.colorbar(p2, ax=ax[1])
cb.set_label('T dust (k)')
ax[1].set_title('Y axis cut')
ax[1].axis('off')
ax[1].set_aspect(1)
p3 = ax[2].contourf(X,Y,avg_Tdust[:,:,nz//2],levels=levels, norm=LogNorm())
cb=f.colorbar(p3, ax=ax[2])
cb.set_label('T dust (k)')
ax[2].set_title('Z axis cut')
ax[2].axis('off')
ax[2].set_aspect(1)
plt.tight_layout()
plt.savefig('figs/AvgDust_T_spec_'+str(species)+'.png', dpi=300)






























