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
from xnemogcm import open_domain_cfg, open_nemo
import cartopy.crs as ccrs
import cmocean
from utils import plot_bathy, plot_bot_plume

# ==============================================================================
# Input parameters

# 1. INPUT FILES

old_bathy = '/data/users/frsy/localdata/UKGO/GO8/bathymetry/old_ORCA025_bathy/domcfg_eORCA025_v2_BATHY_FillZero.nc'
new_bathy = '/data/users/frsy/localdata/UKGO/GO8/bathymetry/new_ORCA025_bathy/BATHY_FROM_GEBCO2021_OCT2021/domain_cfg_eORCA025_noclosea_from_GEBCO2021_S21TT_BATHY.nc'

# 3. PLOT
lon0 = -85.
lon1 = -55.
lat0 =  20.
lat1 =  50.
proj = ccrs.Mercator() #ccrs.Robinson()

# ==============================================================================

# Loading domain geometry
ds_old = xr.open_dataset(old_bathy)
ds_new = xr.open_dataset(new_bathy) 

# Extracting only the part of the domain we need

ds_old = ds_old.isel(x=slice(730,1010),y=slice(760,1000))
ds_new = ds_new.isel(x=slice(730,1010),y=slice(760,1000))

bathy_old = xr.where(ds_old.tmask[0,:,:]==1,ds_old.Bathymetry,np.nan)
bathy_new = xr.where(ds_new.tmask[0,:,:]==1,ds_new.Bathymetry,np.nan)

diff = bathy_old - bathy_new

# Plotting DIFFERENCES ----------------------------------------------------------

fig_name = "old-new_bathymetry.png"
fig_path = "./"
lon = ds_old.nav_lon
lat = ds_old.nav_lat
var = diff
colmap = 'RdBu_r'
vmin = -400.
vmax =  400.
cbar_extend = "both"
cbar_label = "Depth [m]"
cbar_hor = 'horizontal'
map_lims = [lon0, lon1, lat0, lat1]

plot_bathy(fig_name, fig_path, lon, lat, var, proj, 
           colmap, vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims)
