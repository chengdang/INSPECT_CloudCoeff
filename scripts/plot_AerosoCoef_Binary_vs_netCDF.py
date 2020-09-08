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




# Read aerosol coefficint in netCDF format from the netCDF file directly
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




# Read aerosol coefficient obtain from CRTM IO output

# binary I/O
# This folder contains the following output acquied from CRTM Binary LUT:
# Asymmery_factor.txt
# Legendre_terms.txt
# Wavelength.txt
# Extinction_Coefficients.txt
# Radii.txt
# General_Info.txt
# SingleScatAlbedo.txt
crtm_bnfile = '../build/output_aerosol/Binary/'
[bn_iAero, bn_iWvl, bn_iRad] = np.loadtxt(crtm_bnfile+'netCDF_information.txt'). astype(int)
crtm_bn_Reff = np.loadtxt(crtm_bnfile+'Radii.txt')
crtm_bn_pcoef = np.loadtxt(crtm_bnfile+'Legendre_terms.txt')
crtm_bn_Wavelength = np.loadtxt(crtm_bnfile+'Wavelength.txt')
crtm_bn_g  = np.reshape(np.loadtxt(crtm_bnfile+'Asymmery_factor.txt'), (36,61))
crtm_bn_ke = np.reshape(np.loadtxt(crtm_bnfile+'Extinction_Coefficients.txt'), (36,61))
crtm_bn_w  = np.reshape(np.loadtxt(crtm_bnfile+'SingleScatAlbedo.txt'), (36,61))

# netCDF I/O
# This folder contains the following output acquied from CRTM Binary LUT:
# Asymmery_factor.txt
# Legendre_terms.txt
# Wavelength.txt
# Extinction_Coefficients.txt
# Radii.txt
# General_Info.txt
# SingleScatAlbedo.txt
crtm_ncfile = '../build/output_aerosol/NetCDF/'
[nc_iAero, nc_iWvl, nc_iRad] = np.loadtxt(crtm_ncfile+'netCDF_information.txt'). astype(int)
crtm_nc_Reff = np.loadtxt(crtm_bnfile+'Radii.txt')
crtm_nc_pcoef = np.loadtxt(crtm_bnfile+'Legendre_terms.txt')
crtm_nc_Wavelength = np.loadtxt(crtm_bnfile+'Wavelength.txt')
crtm_nc_g  = np.reshape(np.loadtxt(crtm_bnfile+'Asymmery_factor.txt'), (36,61))
crtm_nc_ke = np.reshape(np.loadtxt(crtm_bnfile+'Extinction_Coefficients.txt'), (36,61))
crtm_nc_w  = np.reshape(np.loadtxt(crtm_bnfile+'SingleScatAlbedo.txt'), (36,61))


# Plot
fig = plt.figure(figsize = (12,10))
clines = plt.cm.PiYG(np.linspace(0,1,n_Radii))
# radii
ax = fig.add_subplot(3,2,1)
ax.plot(np.linspace(1, len(crtm_bn_Reff), len(crtm_bn_Reff)), crtm_bn_Reff, '-c', label = 'CRTM_Binary')
ax.plot(np.linspace(1, len(crtm_nc_Reff), len(crtm_nc_Reff)), crtm_nc_Reff, 'om', label = 'CRTM_NetCDF')
ax.plot(np.linspace(1, n_Radii, n_Radii), nc_Reff[bn_iAero-1,:], '.k', label = 'Python NetCDF')
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(labelsize=12)
ax.set_xlabel('#', fontsize = 14)
ax.set_ylabel('radius ($\mathrm{\mu m}$)', fontsize = 14)
plt.legend(loc=2, prop={'size': 12},ncol=1)


# wavelength
ax = fig.add_subplot(3,2,2)
ax.plot(np.linspace(1, len(crtm_bn_Wavelength), len(crtm_bn_Wavelength)), crtm_bn_Wavelength, '-c')
ax.plot(np.linspace(1, len(crtm_nc_Wavelength), len(crtm_nc_Wavelength)), crtm_nc_Wavelength, 'om')
ax.plot(np.linspace(1, n_Wavelengths, n_Wavelengths), nc_Wavelength, '.k')
ax.set_xlabel('#', fontsize = 14)
ax.set_ylabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)

# single-scatteirng albedo
ax = fig.add_subplot(3,2,3)
# for ir in range(n_Radii):
#     ax.plot(crtm_bn_Wavelength,crtm_bn_w[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_w[bn_iAero-1,ir,:],'.k')

ax.plot(crtm_bn_Wavelength,crtm_bn_w[bn_iRad-1,:], '-c')
ax.plot(crtm_nc_Wavelength,crtm_nc_w[nc_iRad-1,:], 'om')
ax.plot(nc_Wavelength,nc_w[bn_iAero-1,bn_iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('single-scattering albedo', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# extiction coefficient
ax = fig.add_subplot(3,2,4)
# for ir in range(n_Radii):
#     ax.plot(crtm_bn_Wavelength,crtm_bn_ke[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_ke[bn_iAero-1,ir,:],'.k')
ax.plot(crtm_bn_Wavelength,crtm_bn_ke[bn_iRad-1,:], '-c')
ax.plot(crtm_nc_Wavelength,crtm_nc_ke[nc_iRad-1,:], 'om')
ax.plot(nc_Wavelength,nc_ke[bn_iAero-1,bn_iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('extinction coefficient ($\mathrm{m^2 kg^-1}$)', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# Aysmmetry factor
ax = fig.add_subplot(3,2,5)
# for ir in range(n_Radii):
#     ax.plot(crtm_bn_Wavelength,crtm_bn_g[ir,:], '-', color = clines[ir])
#     ax.plot(nc_Wavelength,nc_g[bn_iAero-1,ir,:],'.k')
ax.plot(crtm_bn_Wavelength,crtm_bn_g[bn_iRad-1,:], '-c')
ax.plot(crtm_nc_Wavelength,crtm_nc_g[nc_iRad-1,:], 'om')
ax.plot(nc_Wavelength,nc_g[bn_iAero-1,bn_iRad-1,:],'.k')
ax.set_xlabel('wavelength ($\mathrm{\mu m}$)', fontsize = 14)
ax.set_ylabel('asymmetry factor', fontsize = 14)
ax.set_xlim(0.16, 50)
ax.set_xscale('log')

# Phase element
ax = fig.add_subplot(3,2,6)
strwvl = str(round(crtm_bn_Wavelength[bn_iWvl-1],2))
strrad = str(round(crtm_bn_Reff[bn_iRad-1],2))
ax.plot(np.linspace(1, len(crtm_bn_pcoef), len(crtm_bn_pcoef)), crtm_bn_pcoef, '-+c')
ax.plot(np.linspace(1, len(crtm_nc_pcoef), len(crtm_nc_pcoef)), crtm_nc_pcoef, 'om')
nctmp = nc_pcoeff[0, :, bn_iAero-1, bn_iRad-1, bn_iWvl-1]
ax.plot(np.linspace(1, len(nctmp), len(nctmp)), nctmp, '.k')
ax.set_xlabel('Legendre terms' +
              ' ($\mathrm{\lambda}$ = ' + strwvl + '$\mathrm{\mu m}$' +
              ', r = ' + strrad + '$\mathrm{\mu m}$)',
              fontsize = 14)

st = plt.suptitle(Aerosol_Type_String[bn_iAero-1] + ', radius = ' + strrad + '$\mathrm{\mu m}$' , fontsize=22)
st.set_position([.5, 0.95])


plt.show()
