#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 13:58:18 2025

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt

# dirty ice from https://ui.adsabs.harvard.edu/abs/1993A%26A...279..577P/abstract

astrosilfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/sil-dlee.nk'
# astrosilfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/astroSilD2003.nk'
h2ofile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_crystalline.nk'
h2ofile08 = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_Warren08.nk'
dirtyfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/dirty-ice.nk'


astrosildata = np.loadtxt(astrosilfile,max_rows=837,skiprows=2)
h2Odata = np.loadtxt(h2ofile,max_rows=2336,skiprows=2)
h2O08data = np.loadtxt(h2ofile08,max_rows=486,skiprows=2)
dirtydata = np.loadtxt(dirtyfile,max_rows=90,skiprows=2)

plt.figure()
plt.xscale('log')
# plt.yscale('log')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(h2Odata[:,0], h2Odata[:,1],'C1',label='h2O ice n')
plt.plot(h2Odata[:,0], h2Odata[:,2],'C1',label='h2O ice k')

plt.plot(h2O08data[:,0], h2O08data[:,1],'C2--',label='h2O ice new n')
plt.plot(h2O08data[:,0], h2O08data[:,2],'C2--',label='h2O ice new k')


plt.plot(dirtydata[:,0], dirtydata[:,1],'C3--',label='dirty ice n')
plt.plot(dirtydata[:,0], dirtydata[:,2],'C3--',label='dirty ice k')

plt.legend()

newlambda = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,0], h2O08data[:,0]))
new_n = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,1], h2O08data[:,1]))
new_k = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,2], h2O08data[:,2]))

# plt.plot(newlambda, new_n,'C2',label='n')
# plt.plot(newlambda, new_k,'C2',label='k')


header = "nk\nh20_ice.nk 273 0.917 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/H2O_ice.nk'
np.savetxt(output_file, np.vstack((newlambda, new_n, new_k)).T, header=header, comments='', fmt='%.6E')




















