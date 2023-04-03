#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as feature
from utils_sec import *

# ==============================================================================
# Input parameters

# 1. INPUT FILES

config = 'orca12'
maxdep = 500. #5800. 

if config == 'orca025':
   bathy_file = '/data/users/dbruciaf/GS/orca025/GO8_domain/bathy_eORCA025_no_closea-20220404.nc'
   coord_file = '/data/users/dbruciaf/GS/orca025/GO8_domain/coordinates.nc'
else:
   bathy_file = '/data/users/dbruciaf/GS/orca12/GO8_domain/bathy_eORCA12_no_closea.nc'
   coord_file = '/data/users/dbruciaf/GS/orca12/GO8_domain/coordinates.nc'

sec_lon1 = [-81.11, -70.38]
sec_lat1 = [ 28.15,  27.80]
sec_lon2 = [-81.85, -73.46]
sec_lat2 = [ 29.92,  30.56]
sec_lon3 = [-81.33, -73.28]
sec_lat3 = [ 32.06,  31.65]
sec_lon4 = [-82.56, -66.82]
sec_lat4 = [ 32.82,  26.69]
sec_lon5 = [-80.01, -69.55]
sec_lat5 = [ 33.14,  32.44]
sec_lon6 = [-78.02, -67.90]
sec_lat6 = [ 34.62,  32.25]
# Cape Hatteras
sec_lon7 = [-75.91, -67.38]
sec_lat7 = [ 35.69,  34.33]
#
sec_lon8 = [-75.58, -67.38]
sec_lat8 = [ 38.56,  35.42]
# Ezer 2016 sec
sec_lon9 = [-77., -60.]
sec_lat9 = [ 35.,  35.]
# North Cape Hatteras
sec_lon10 = [-78., -60.]
sec_lat10 = [ 37.,  37.]
# 
sec_lon11 = [-70.11, -65.54] 
sec_lat11 = [ 44.40,  38.83]
#
sec_lon12 = [-58.99, -54.93, -47.96]
sec_lat12 = [ 51.75,  48.18,  39.15]



sec_i = [sec_lon1, sec_lon2, sec_lon3, sec_lon4, sec_lon5, sec_lon6, sec_lon7, sec_lon8, sec_lon9, sec_lon10, sec_lon11, sec_lon12]
sec_j = [sec_lat1, sec_lat2, sec_lat3, sec_lat4, sec_lat5, sec_lat6, sec_lat7, sec_lat8, sec_lat9, sec_lat10, sec_lat11, sec_lat12]

proj = ccrs.Mercator()
transform = ccrs.PlateCarree()
# ==============================================================================

# Loading domain geometry
ds_bat  = xr.open_dataset(bathy_file)
ds_cor  = xr.open_dataset(coord_file)

# Extracting only the part of the domain we need
if config == 'orca025':
   ds_bat = ds_bat.isel(x=slice(730,1100), y=slice(760,1100))
   ds_cor = ds_cor.isel(x=slice(730,1100), y=slice(760,1100))
else:
   ds_bat = ds_bat.isel(x=slice(2230,3100), y=slice(2065,3300))
   ds_cor = ds_cor.isel(x=slice(2230,3100), y=slice(2065,3300))

# Plotting BATHYMETRY ----------------------------------------------------------

bathy = ds_bat["Bathymetry"].squeeze() #.isel(x_c=slice(1, None), y_c=slice(1, None))
lon = ds_cor["glamf"].squeeze()
lat = ds_cor["gphif"].squeeze()

for s in range(len(sec_i)):
    I = sec_i[s]
    J = sec_j[s]
    secI = []
    secJ = []
    for p in range(len(I)):
        j, i = get_ij_from_lon_lat(I[p], J[p], lon.data, lat.data)
        secI.append(i)
        secJ.append(j)
    SEC_J, SEC_I = get_poly_line_ij(np.asarray(secI), np.asarray(secJ))
    DEP = []
    LON = []
    LAT = []
    for z in range(len(SEC_I)):
        DEP.append(bathy.data[SEC_J[z], SEC_I[z]])
        LON.append(lon.data[SEC_J[z], SEC_I[z]])
        LAT.append(lat.data[SEC_J[z], SEC_I[z]])
    DEP = np.asarray(DEP) 
    LON = np.asarray(LON)
    LAT = np.asarray(LAT)

    fig = plt.figure(figsize = (34.,17.5), dpi=100)
    ax = fig.gca()
    zd = plt.plot(DEP, color='red')
    plt.setp(zd, 'linewidth', 4., zorder=20)
    ax.set_ylim([0., maxdep])
    ax.invert_yaxis()    

    #ax.set_xlabel(xlabel, fontsize='50')
    ax.set_ylabel('Depth [m]', fontsize='50')
    ax.tick_params(axis='x', labelsize=50)
    ax.tick_params(axis='y', labelsize=50)

    pos1 = ax.get_position()
    lon0 = -85.
    lon1 = -30.
    lat0 =  20.
    lat1 =  65.
 
    map_lims = [lon0, lon1, lat0, lat1]
    a = plt.axes([pos1.x0-0.002, pos1.y0+0.01, .3, .3], projection=proj)
    a.coastlines()
    a.add_feature(feature.LAND, color='gray',edgecolor='gray',zorder=1)
    MAP = plt.plot(LON, LAT, c="red", transform=transform)
    plt.setp(MAP, 'linewidth', 2.5)
    a.set_extent(map_lims)

    fig_name = 'sec_lon_'+str(LON[0])+'_lat_'+str(LAT[0])+'_'+str(maxdep)+'.png'
    fig_path = './'
    print(f"Saving {fig_path}", end=": ")
    plt.savefig(fig_path+fig_name, bbox_inches="tight")
    print("done")
    plt.close()
 
