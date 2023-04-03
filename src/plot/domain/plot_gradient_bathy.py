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
from utils_grad import plot_bathy, plot_bot_plume
from scipy.ndimage import gaussian_filter

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
ds_cor  = xr.open_dataset(coord_file).squeeze()

# Extracting only the part of the domain we need

ds_bat = ds_bat.isel(x=slice(700,1500), y=slice(760,1300))
ds_cor = ds_cor.isel(x=slice(700,1500), y=slice(760,1300))

H_T = ds_bat["Bathymetry"]
H_T = H_T.where(H_T>1000.,0.)
H_T = gaussian_filter(H_T, sigma=2.)
#e1t = ds_cor.e1t
#e2t = ds_cor.e2t
e1u = ds_cor.e1u
e2v = ds_cor.e2v

# Computing gradient components
dHdx = np.diff(H_T, axis=1, append=np.nan) / e1u # gradient is now on U points
dHdx = dHdx.rolling({'x':2}).mean().fillna(0.) # back to T points

dHdy = np.diff(H_T, axis=0, append=np.nan) / e2v # gradient is now on V points
dHdy = dHdy.rolling({'y':2}).mean().fillna(0.) # back to T points

grad_H = np.sqrt(dHdx**2 + dHdy**2)
grad_H = grad_H.where(ds_bat["Bathymetry"]>900.)
print(np.nanmax(grad_H))

# Plotting gradient ----------------------------------------------------------


fig_name = "bathymetry_gradient.png"
fig_path = "./"
lon = ds_cor["glamf"].squeeze()
lat = ds_cor["gphif"].squeeze()
var1 = grad_H
var2 = ds_bat["Bathymetry"] 
colmap = cmocean.cm.deep
vmin = 0.
vmax = 6000.
cbar_extend = "max"
cbar_label = "Depth [m]"
cbar_hor = 'horizontal'
map_lims = [lon0, lon1, lat0, lat1]

plot_bathy(fig_name, fig_path, lon, lat, var1, var2, proj, 
           colmap, vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims)

 
