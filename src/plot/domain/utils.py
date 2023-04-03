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
import matplotlib.gridspec as gridspec
import xarray as xr
from xnemogcm import open_domain_cfg, open_nemo
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import cmocean

def plot_bathy(fig_name, fig_path, lon, lat, var, proj, colmap, vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims):

    fig = plt.figure(figsize=(20,20), dpi=100)
    spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

    ax = fig.add_subplot(spec[:1], projection=proj)
    ax.coastlines()
    #ax.gridlines()

    # Drawing settings
    cmap = colmap #cmocean.cm.deep
    transform = ccrs.PlateCarree()
    pcol_kwargs = dict(cmap=cmap, vmin=vmin, vmax=vmax, transform=transform)
    cbar_kwargs = dict(extend=cbar_extend, orientation=cbar_hor)
    #plot_kwargs = dict(color="r", transform=ccrs.PlateCarree())
    land_col = ".15"

    # Grid settings
    gl_kwargs = dict()
    gl = ax.gridlines(**gl_kwargs)
    gl.xlines = False
    gl.ylines = False
    gl.top_labels = True
    gl.right_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 30, 'color': 'k'}
    gl.ylabel_style = {'size': 30, 'color': 'k'}

    # Plotting
    #LON = [-81., -83.5, -74. , -65.5, -65.5, -81.]
    #LAT = [ 26.,  35. ,  43.6,  38.0,  26. ,  26.]
    #LON = [-81., -83.5, -71.5 , -70.5, -53.0, -41.5,-41.5,-81.0]
    #LAT = [ 26.,  35. ,  45.5,   56.0,  60.0 , 47.0, 26.0, 26.0]

    #LON = [-81., -83.5, -71.5 , -70.5, -48.0, -40.0,-53.0,-81.0]
    #LAT = [ 26.,  35. ,  45.5,   50.0,  56.0,  47.0, 26.0, 26.0]
    
    ##LON = [-81., -83.5, -71.5 , -71.0, -46.0, -46.0, -38.0, -38.0,-51.0,-81.0]
    ##LAT = [ 26.,  35. ,  45.5,   67.0,  66.0,  61.5,  57.0,  47.0, 26.0, 26.0]
    LON = [-81., -83.5, -71.5 , -71.0, -46.0, -46.0, -39.0, -41.2, -39.0,-51.0,-81.0]
    LAT = [ 26.,  35. ,  45.5,   67.0,  66.0,  61.5,  58.0,  50.0, 47.0, 26.0, 26.0]
    var = var.where(np.isfinite(var), -1)
    bathy = var 
    if 'x_c' in var.dims:
       var = var.isel(x_c=slice(1, None), y_c=slice(1, None))
    else:
       var = var.isel(x=slice(1, None), y=slice(1, None))
    pcol = ax.pcolormesh(lon, lat, var, **pcol_kwargs)
    pcol.cmap.set_under(land_col)
    bcon = ax.contour(lon, lat, bathy, levels=[300.], colors='magenta', linewidths=2, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[400.], colors='black', linewidths=2, transform=transform)
    #bcon = ax.contour(lon, lat, bathy, levels=[3000.], colors='magenta', linewidths=2, transform=transform)
    #bcon = ax.contour(lon, lat, bathy, levels=[3500.], colors='black', linewidths=2, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[4350.], colors='gold', linewidths=4, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[4000.], colors='cyan', linewidths=2, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[4500.], colors='limegreen', linewidths=2, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[5100.], colors='red', linewidths=2, transform=transform)
    bcon = ax.contour(lon, lat, bathy, levels=[0.002], colors='magenta', linewidths=2, transform=transform)
    ax.plot(LON, LAT, 'o-',color='red', linewidth=5, transform=transform)

    cb = plt.colorbar(pcol, **cbar_kwargs)
    cb.set_label(label=cbar_label,size=40)
    cb.ax.tick_params(labelsize=30)
    ax.set_extent(map_lims)
    print(f"Saving {fig_path}", end=": ")
    plt.savefig(fig_path+fig_name, bbox_inches="tight")
    print("done")
    plt.close()


def plot_bot_plume(fig_name, fig_path, lon, lat, var, proj, colmap, vmin, vmax, cbar_extend, cbar_label, cbar_hor, map_lims, bathy=None, cn_lev=None):

    fig = plt.figure(figsize=(20,20), dpi=100)
    spec = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)

    ax = fig.add_subplot(spec[:1], projection=proj)
    #ax.coastlines()
    #ax.gridlines()

    # Drawing settings
    cmap = colmap #cmocean.cm.deep
    transform = ccrs.PlateCarree()
    pcol_kwargs = dict(cmap=cmap, vmin=vmin, vmax=vmax, transform=transform)
    cbar_kwargs = dict(extend=cbar_extend, orientation=cbar_hor)
    #plot_kwargs = dict(color="r", transform=ccrs.PlateCarree())
    land_col = ".15"

    # Grid settings
    gl_kwargs = dict()
    gl = ax.gridlines(**gl_kwargs)
    gl.xlines = False
    gl.ylines = False
    gl.top_labels = True
    gl.right_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 30, 'color': 'k'}
    gl.ylabel_style = {'size': 30, 'color': 'k'}

    # Plotting
    var = var.where(var != 0, -1)
    land = -1. + bathy.where(bathy == 0)
    pcol = ax.pcolormesh(lon, lat, var, alpha=0.8, **pcol_kwargs)
    if bathy is not None:
       bcon = ax.contour(lon, lat, bathy, levels=cn_lev, colors='k', transform=transform)
       bcon = ax.contour(lon, lat, bathy, levels=[2300.], colors='r', linewidths=5.0, transform=transform)
       bcol = ax.pcolormesh(lon, lat, land, **pcol_kwargs)
       bcol.cmap.set_under(land_col)
    cb = plt.colorbar(pcol, **cbar_kwargs)
    cb.set_label(label=cbar_label,size=40)
    cb.ax.tick_params(labelsize=30) 
    ax.set_extent(map_lims)
    print(f"Saving {fig_path}", end=": ")
    plt.savefig(fig_path+fig_name, bbox_inches="tight")
    print("done")
    plt.close()
