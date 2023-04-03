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
import cmocean
from utils import plot_bathy, plot_bot_plume

# ==============================================================================
# Input parameters

# 1. INPUT FILES

bathy_file = '/data/users/dbruciaf/GS/orca025/GO8_domain/bathy_eORCA025_no_closea-20220404.nc'
coord_file = '/data/users/dbruciaf/GS/orca025/GO8_domain/coordinates.nc'

# 3. PLOT
lon0 = -95.
lon1 = -10.
lat0 =  20.
lat1 =  70.
proj = ccrs.Mercator() #ccrs.Robinson()

# ==============================================================================

# Loading domain geometry
ds_bat  = xr.open_dataset(bathy_file)
ds_cor  = xr.open_dataset(coord_file)

# Extracting only the part of the domain we need

ds_bat = ds_bat.isel(x=slice(700,1500), y=slice(760,1300))
ds_cor = ds_cor.isel(x=slice(700,1500), y=slice(760,1300))

# Plotting BATHYMETRY ----------------------------------------------------------

bathy = ds_bat["Bathymetry"] #.isel(x_c=slice(1, None), y_c=slice(1, None))


fig_name = "bathymetry.png"
fig_path = "./"
lon = ds_cor["glamf"].squeeze()
lat = ds_cor["gphif"].squeeze()
var = bathy
colmap = cmocean.cm.deep
vmin = 0.
vmax = 6000.
cbar_extend = "max"
cbar_label = "Bathymetry [m]"
cbar_hor = 'horizontal'
map_lims = [lon0, lon1, lat0, lat1]

plot_bathy(fig_name, fig_path, lon, lat, var, proj, 
           colmap, vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims)

 
