#!/usr/bin/env python
# coding: utf-8


# Process Aerosol SSP
# Binary vs. netCDF
# Cheng Dang, Aug 1, 2020


from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

def readnc(infile, varname):
    nc_fid = Dataset(infile,'r')
    out =  np.array(nc_fid.variables[varname][:])
    return out

Aerosol_Type_String = ['DUST','SEASALT_SSAM','SEASALT_SSCM1','SEASALT_SSCM2',
                       'SEASALT_SSCM3','ORGANIC_CARBON','BLACK_CARBON','SULFATE']




# Read aerosol coefficint in netCDF format
ncfile = '../fix/AerosolCoeff.nc'
n_Wavelengths = 61 ;
n_Radii = 36 ;
n_Types = 8 ;
n_RH = 36 ;
n_Legendre_Terms = 38 ;
n_Phase_Elem = 1;
nc_Aerosol_Type = readnc(ncfile,'Aerosol_Type')
nc_Aerosol_Type_Name = readnc(ncfile,'Aerosol_Type_Name')
nc_Wavelength = readnc(ncfile,'Wavelength')
nc_Reff = readnc(ncfile,'Reff')
nc_RH = readnc(ncfile,'RH')
nc_ke = readnc(ncfile,'ke')
nc_w  = readnc(ncfile,'w')
nc_g  = readnc(ncfile,'g')
nc_pcoeff = readnc(ncfile,'pcoeff')




# Read aerosol coefficient obtain from binary I/O
# This folder contains the following output acquied from CRTM Binary LUT:
# Asymmery_factor.txt
# Legendre_terms.txt
# Wavelength.txt
# Extinction_Coefficients.txt
# Radii.txt
# General_Info.txt
# SingleScatAlbedo.txt

bnfile = '../build/output_aerosol/'
[iAero, iWvl, iRad] = np.loadtxt(bnfile+'netCDF_information.txt'). astype(int)
bn_Reff = np.loadtxt(bnfile+'Radii.txt')
bn_pcoef = np.loadtxt(bnfile+'Legendre_terms.txt')
bn_Wavelength = np.loadtxt(bnfile+'Wavelength.txt')
bn_g  = np.reshape(np.loadtxt(bnfile+'Asymmery_factor.txt'), (36,61))
bn_ke = np.reshape(np.loadtxt(bnfile+'Extinction_Coefficients.txt'), (36,61))
bn_w  = np.reshape(np.loadtxt(bnfile+'SingleScatAlbedo.txt'), (36,61))



# Plot
fig = plt.figure(figsize = (16,16))
clines = plt.cm.PiYG(np.linspace(0,1,n_Radii))
# radii
ax = fig.add_subplot(3,2,1)
ax.plot(np.linspace(1, len(bn_Reff), len(bn_Reff)), bn_Reff, '-c', label = 'Binary')
ax.plot(np.linspace(1, n_Radii, n_Radii), nc_Reff[iAero-1,:], '.k', label = 'NetCDF')
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(labelsize=12)
ax.set_xlabel('#', fontsize = 14)
ax.set_ylabel('radius ($\mathrm{\mu m}$)', fontsize = 14)
plt.legend(loc=2, prop={'size': 12},ncol=1)


# wavelength
ax = fig.add_subplot(3,2,2)
ax.plot(np.linspace(1, len(bn_Wavelength), len(bn_Wavelength)), bn_Wavelength, '-c')
ax.plot(np.linspace(1, n_Wavelengths, n_Wavelengths), nc_Wavelength, '.k')
ax.set_xlabel('#', fontsize = 14)
ax.set_ylabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)

# single-scatteirng albedo
ax = fig.add_subplot(3,2,3)
# for ir in range(n_Radii):
#     ax.plot(bn_Wavelength,bn_w[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_w[iAero-1,ir,:],'.k')

ax.plot(bn_Wavelength,bn_w[iRad-1,:], '-c')
ax.plot(nc_Wavelength,nc_w[iAero-1,iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('single-scattering albedo', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# extiction coefficient
ax = fig.add_subplot(3,2,4)
# for ir in range(n_Radii):
#     ax.plot(bn_Wavelength,bn_ke[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_ke[iAero-1,ir,:],'.k')
ax.plot(bn_Wavelength,bn_ke[iRad-1,:], '-c')
ax.plot(nc_Wavelength,nc_ke[iAero-1,iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('extinction coefficient ($\mathrm{m^2 kg^-1}$)', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# Aysmmetry factor
ax = fig.add_subplot(3,2,5)
# for ir in range(n_Radii):
#     ax.plot(bn_Wavelength,bn_g[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_g[iAero-1,ir,:],'.k')
ax.plot(bn_Wavelength,bn_g[iRad-1,:], '-c')
ax.plot(nc_Wavelength,nc_g[iAero-1,iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('asymmetry factor', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# Phase element
ax = fig.add_subplot(3,2,6)
strwvl = str(round(bn_Wavelength[iWvl-1],2))
strrad = str(round(bn_Reff[iRad-1],2))
ax.plot(np.linspace(1, len(bn_pcoef), len(bn_pcoef)), bn_pcoef, '-+c')
nctmp = nc_pcoeff[0, :, iAero-1, iRad-1, iWvl-1]
ax.plot(np.linspace(1, len(nctmp), len(nctmp)), nctmp, '.k')
ax.set_xlabel('Legendre terms' +
              ' ($\mathrm{\lambda}$ = ' + strwvl + '$\mathrm{\mu m}$' +
              ', r = ' + strrad + '$\mathrm{\mu m}$)',
              fontsize = 14)

st = plt.suptitle(Aerosol_Type_String[iAero-1] + ', radius = ' + strrad + '$\mathrm{\mu m}$' , fontsize=22)
st.set_position([.5, 0.95])

plt.show()


