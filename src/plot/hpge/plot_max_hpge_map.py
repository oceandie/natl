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
from xnemogcm import open_domain_cfg
import cartopy.crs as ccrs
import cmocean
from utils import plot_hpge

# ==============================================================================
# Input parameters

# 1. INPUT FILES

vcoord = 'MEs'
DOMCFG_file = '/data/users/dbruciaf/GS/orca025/MEs_GO8/4env_MEs_000_015_020/domain_cfg_000-015-020_v1.nc'
HPGE_dir = '/data/users/dbruciaf/GS/orca025/MEs_GO8/4env_MEs_000_015_020/'
HPGE_file = 'maximum_hpge_000-015-020_v1_tot.nc'

# 3. PLOT
lon0 = -90.
lon1 = -35.
lat0 =  20.
lat1 =  60.
proj = ccrs.Mercator() #ccrs.Robinson()

# ==============================================================================

# Loading domain geometry
ds_dom  = open_domain_cfg(files=[DOMCFG_file])
for i in ['bathymetry','bathy_meter']:
    for dim in ['x','y']:
        ds_dom[i] = ds_dom[i].rename({dim: dim+"_c"})

ds_hpge  = xr.open_dataset(HPGE_dir + HPGE_file)

# Extracting only the part of the domain we need

ds_dom = ds_dom.isel(x_c=slice(730,1100),x_f=slice(730,1100),
                     y_c=slice(760,1200),y_f=slice(760,1200))
ds_hpge = ds_hpge.isel(x=slice(730,1100),y=slice(760,1200))

# Plotting BATHYMETRY ----------------------------------------------------------

#bathy = ds_dom["bathymetry"]
bathy = ds_dom["bathymetry"]#.isel(x_c=slice(1, None), y_c=slice(1, None))
varss = list(ds_hpge.keys())

for env in range(len(varss)):

    fig_name = varss[env] + '_' + vcoord + '.png'
    fig_path = "./"
    lon = ds_dom["glamf"]
    lat = ds_dom["gphif"]
    var = ds_hpge[varss[env]] 
    colmap = 'hot' #cmocean.cm.ice
    vmin = 0.0
    vmax = 0.05
    cbar_extend = 'max' #"max"
    cbar_label = "HPG errors [$m\;s^{-1}$]"
    cbar_hor = 'horizontal'
    map_lims = [lon0, lon1, lat0, lat1]
    #cn_lev = [0., 250., 500., 750., 1000., 1250., 1500, 1750., 2000., 2500., 3000., 3500.] 
    cn_lev = [500.]

    plot_hpge(fig_name, fig_path, lon, lat, var, proj, colmap, 
              vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims, bathy, cn_lev)

 
