import numpy as np
import matplotlib.pyplot as plt
import sys, os

# add mocassin tools module
sys.path.insert(0, 'mocassin_fit')
from mocassin_tools import *

from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.visualization import make_lupton_rgb
import scipy.ndimage as ndimage
from scipy.signal import savgol_filter, medfilt
import warnings
warnings.filterwarnings('ignore')

plt.close('all')

############################################################
# define model directory

model_dir = '/home/hmonteiro/software/mocassin_HM_2026/benchmarks/test_Problem/'

mod_id = model_dir.split(sep='/')[-2]

############################################################
# Check the existence of important directories for output

if not os.path.exists(model_dir+'figs/'):
    os.makedirs(model_dir+'figs/')
    print(f"Diretório '{model_dir+'figs/'}' criado com sucesso!")
else:
    print(f"O diretório '{model_dir+'figs/'}' já existe.")
    
if not os.path.exists(model_dir+'maps/'):
    os.makedirs(model_dir+'maps/')
    print(f"Diretório '{model_dir+'maps/'}' criado com sucesso!")
else:
    print(f"O diretório '{model_dir+'maps/'}' já existe.")
    
############################################################
# Read input file from mocassin

print ('Reading in input parameters from input.in ...')
moc_in = open(model_dir+"input/input.in", "r")
lines = moc_in.readlines()
moc_pars = []
for line in lines:
    if (line.split() != []):
        moc_pars.append(line.split())
moc_in.close()    

if ['symmetricXYZ'] in moc_pars: 
    do_x8 = True
else:
    do_x8 = False

# define inclination angle
incl_angle = 40.

############################################################
# Reading Grid Structure
        
print ('Reading Grid Structure ...')
grid0 = open(model_dir+"output/grid0.out", "r")
lines = grid0.readlines()

nx,ny,nz, aux = [int(x) for x in (lines[1].split())[0:4]]
ncell = nx*ny*nz

aux,Rout = [float(x) for x in (lines[1].split())[4:]]
X = np.array([float(x) for x in lines[2:2+nx]])
Y = np.array([float(x) for x in lines[2+nx:2+nx+ny]])
Z = np.array([float(x) for x in lines[2+nx+ny:2+nx+ny+nz]])

XX, YY, ZZ = np.meshgrid(X,Y,Z, indexing='ij')

aux = np.array([[int(y) for y in x.split()] for x in lines[2+nx+ny+nz:]])
converged = np.reshape(aux[:,1],(nx,ny,nz))

plt.figure()
plt.imshow(converged[:,:,int(nx/2)])
plt.title('Cell convergente')

if (nx == ny == nz or Rout == 0):
    Rout = np.abs(Z).max() 

############################################################
# reading  electron temperatures, electron densities and hydrogen densities

print ('Reading Ne and Te Structure ...')
grid1 = np.loadtxt(model_dir+'output/grid1.out')

Te = np.reshape(grid1[:,0],(nx,ny,nz))
Ne = np.reshape(grid1[:,1],(nx,ny,nz))
H0 = np.reshape(grid1[:,2],(nx,ny,nz))

if do_x8:
    Te = cube_rot(make_cube(Te),incl_angle,0,0,order=1)
    Ne = cube_rot(make_cube(Ne),incl_angle,0,0,order=1)
    H0 = cube_rot(make_cube(H0),incl_angle,0,0,order=1)


print('Peak density in H0: ',H0.max())

f,ax = plt.subplots(1,3,figsize=(11,3))
p1 = ax[0].imshow(H0.mean(axis=0), origin='lower')
f.colorbar(p1, ax=ax[0])
ax[0].set_title('H density ($cm^{-3}$)')
ax[0].axis('off')
p2 = ax[1].imshow(Ne.mean(axis=0), origin='lower')
f.colorbar(p2, ax=ax[1])
ax[1].set_title('Ne ($cm^{-3}$)')
ax[1].axis('off')
p3 = ax[2].imshow(Te.mean(axis=0), origin='lower')
f.colorbar(p3, ax=ax[2])
ax[2].set_title('Te (K)')
ax[2].axis('off')
plt.tight_layout()
plt.savefig(model_dir+'figs/projected_H0_Ne_Te_mod.png', dpi=300)

ind_cut = H0.shape[0]//2
    
f,ax = plt.subplots(1,3,figsize=(11,3))
p1 = ax[0].imshow(H0[ind_cut,:,:], origin='lower',norm=LogNorm(vmin=500))
f.colorbar(p1, ax=ax[0])
ax[0].set_title('H density ($cm^{-3}$)')
ax[0].axis('off')
p2 = ax[1].imshow(Ne[ind_cut,:,:], origin='lower',norm=LogNorm(vmin=500))
f.colorbar(p2, ax=ax[1])
#ax[1].contour(H0,colors='w')
ax[1].set_title('Ne ($cm^{-3}$)')
ax[1].axis('off')
p3 = ax[2].imshow(Te[ind_cut,:,:], origin='lower',norm=LogNorm(vmin=8000))
f.colorbar(p3, ax=ax[2])
ax[2].set_title('Te (K)')
ax[2].axis('off')
plt.tight_layout()
plt.savefig(model_dir+'figs/Zcut_H0_Ne_Te_mod.png', dpi=300)

############################################################
# Read input file from mocassin

print ('Reading in lineflux.out file ...')
lineflux = open(model_dir+"output/lineFlux.out", "r")
lines = lineflux.readlines()
line_flux_out = []
for line in lines:
    if (line.split() != []):
        line_flux_out.append(line.split())
lineflux.close()    

hbeta = float(line_flux_out[2][3])

############################################################
# Other data needed
dist = 1000.

mod_hbeta = (hbeta*1.e36)/4./np.pi/(3.086e18*dist)**2
print ('Using distance of d =%5.0f parsecs.'%dist)
print ('Model Hbeta at earth: %6.2e   Lit.: 7.76e-10 '%mod_hbeta)
 
############################################################
# observed size for given distance
arcsec_sz = 3600.*180./np.pi*np.arctan(Rout/(dist*3.086e18))


############################################################
# reading  line emissivities
#%%

# determine arcsec per cell
arcsec_cell = arcsec_sz/(nz/2)

print ('Reading line emissivities ...')
lines = np.loadtxt(model_dir+'input/plot.in', skiprows=1, usecols=2)
emiss = np.loadtxt(model_dir+'output/plot.out')

# generate line maps  and save fits files
plot_lines = [4861., 5008., 6585., 4686., 6718.,1906.]

cont = 1
line_maps = []

for i in range(len(lines)):
#    print ('Doing line %5.0f '%lines[i], i)
    cube = np.reshape(emiss[:,i+1],(nx,ny,nz))
    if do_x8:
        # correct for escaped border photons in symmetricXYZ mode 
        cube[0,:,:] = cube[0,:,:] * 2
        cube[:,0,:] = cube[:,0,:] * 2
        cube[:,:,0] = cube[:,:,0] * 2
        cube[-1,:,:] = cube[-1,:,:] * 2
        cube[:,-1,:] = cube[:,-1,:] * 2
        cube[:,:,-1] = cube[:,:,-1] * 2
        cube = make_cube(cube)
        cube = cube_rot(cube,incl_angle,0,0,order=1)
        


    if(lines[i] == 6563.):
        hacube = cube
    if(lines[i] == 4861.):
        hbcube = cube
    if(lines[i] == 6585.):
        n26583cube = cube
    if(lines[i] == 5756.):
        n25756cube = cube
    if(lines[i] == 9072.):
        s39072cube = cube
    if(lines[i] == 4364.):
        o34363cube = cube
    if(lines[i] == 6302.):
        o16300cube = cube
    line_maps.append(np.nansum(cube,axis=0))

    hdu = fits.PrimaryHDU(np.flip(np.nansum(cube,axis=0),axis=0)*1.e36)
    hdulist = fits.HDUList([hdu])
        
    hdulist[0].header['Lambda'] = str(lines[i])
    hdulist[0].header['CTYPE1'] = 'RA---TAN'
    hdulist[0].header['CTYPE2'] = 'DEC--TAN'
    hdulist[0].header['CDELT1'] = (arcsec_cell*u.arcsec).to(u.degree).value 
    hdulist[0].header['CDELT2'] = (arcsec_cell*u.arcsec).to(u.degree).value
    hdulist[0].header['CUNIT1'] = 'deg'
    hdulist[0].header['CUNIT2'] = 'deg'
    hdulist[0].header['CRVAL1'] =   258.43500          
    hdulist[0].header['CRVAL2'] =   -37.10325 
           
    if do_x8:
        hdulist[0].header['CRPIX1'] = int(nz)+1
        hdulist[0].header['CRPIX2'] = int(ny)+1
    else:
        hdulist[0].header['CRPIX1'] = int(nz/2)+1
        hdulist[0].header['CRPIX2'] = int(ny/2)+1
        
    
    hdulist.writeto(model_dir+'maps/'+mod_id+'_model_'+str(incl_angle)+'_'+str(lines[i])+'.fits',overwrite=True, output_verify='fix')


###############################################################################
# Plot line maps
plot_lines = [4861., 5008., 6585., 4686., 6718.,5877.]

plt.figure()
cont = 1
for line in plot_lines:
    ax = plt.subplot(2,3,cont)
    ind = np.where(lines == line)
    ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(),norm=LogNorm(vmin=1.e-6,vmax=0.99), origin='lower')
    ax.axis('off')
    #ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(), vmin=1.e-8,vmax=0.6)
    ax.set_title(str(line))
    cont += 1
plt.tight_layout()
plt.savefig(model_dir+'figs/model_linemaps-log.png', dpi=300)


plt.figure()
cont = 1
for line in plot_lines:
    ax = plt.subplot(2,3,cont)
    ind = np.where(lines == line)
    ax.axis('off')
    ax.imshow((line_maps[ind[0][0]]+1.0e-16)/line_maps[ind[0][0]].max(), vmin=0.,vmax=0.99, origin='lower')
    ax.set_title(str(line))
    cont += 1
plt.tight_layout()
plt.savefig(model_dir+'figs/model_linemaps.png', dpi=300)

###############################################################################
# Plot Ha/Hb
#%%
indHb = np.where(lines == 4861.)
indHa = np.where(lines == 6563.)

plt.figure()
plt.title(r'H$\alpha$/H$\beta$')
plt.imshow(line_maps[indHa[0][0]]/line_maps[indHb[0][0]])

############################################################
# reading  ionic fractions
#%%
Nelements = 9

H_ionic_frac = np.zeros((2,ncell))
He_ionic_frac = np.zeros((3,ncell))
Heavy_ionic_frac = np.zeros((Nelements-2,7,ncell))

grid2 = open(model_dir+"output/grid2.out", "r")
aux = grid2.readlines()

for i in range(ncell):
    H_ionic_frac[:,i] = aux[i*Nelements].split()
    He_ionic_frac[:,i] = aux[i*Nelements+1].split()
    Heavy_ionic_frac[:,:,i] = np.array([x.split() for x in aux[i*Nelements+2:i*Nelements+2+(Nelements-2)]])

# dealocate memory 
del aux

