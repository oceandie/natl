#!/usr/bin/env pycnd2

import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# MEs domain_cfg,nc
cfg = "/data/users/dbruciaf/OVF/MEs_GO8/env4.v1/domain_cfg.4env.v1.nc"
nc      = nc4.Dataset(cfg,"r")
stf3    = np.array(nc.variables["stiff3D"])[0,:,900:2001,900:1160]
lon     = np.array(nc.variables["nav_lon"])[900:2001,900:1160]
lat     = np.array(nc.variables["nav_lat"])[900:2001,900:1160]
msk     = np.array(nc.variables["top_level"])[0,900:2001,900:1160]
nc.close()

stf = np.ones(shape=lon.shape)*np.nan
lev = np.ones(shape=lon.shape)*np.nan

for j in range(lon.shape[0]):
    for i in range(lon.shape[1]):
        if msk[j,i] == 1:
           stf[j,i] = np.nanmax(stf3[:,j,i])
           lev[j,i] = np.where(stf3[:,j,i]==stf[j,i])[0][0]

LLcrnrlon = -42.50 
LLcrnrlat =  53.00 
URcrnrlon =   4.
URcrnrlat =  73.50

print 'Computing projection ...'
proj = Basemap(projection='merc',resolution='h',
               llcrnrlon=LLcrnrlon,llcrnrlat=LLcrnrlat,
               urcrnrlon=URcrnrlon,urcrnrlat=URcrnrlat)

x2 , y2  = proj(lon, lat) 

fig, ax = plt.subplots(figsize=(50, 25))
proj.drawcoastlines(zorder=4)
proj.fillcontinents('grey', zorder=4)
proj.drawmapboundary()
proj.drawparallels(np.arange(-90,90,5.),labels=[True,False,False,True],fontsize=30)
proj.drawmeridians(np.arange(-180,180,5.),labels=[True,False,False,True],fontsize=30)
pc = proj.pcolormesh(x2, y2, stf, cmap="hot_r",vmin=0., vmax=40., zorder=3)
cbar = proj.colorbar(pc)
cbar.ax.tick_params(labelsize=30)
cbar.set_label('Max stiffness parameter [$m$]',size=40)
out_name ='max_stiff.png'
plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1)
plt.close()

fig, ax = plt.subplots(figsize=(50, 25))
proj.drawcoastlines(zorder=4)
proj.fillcontinents('grey', zorder=4)
proj.drawmapboundary()
proj.drawparallels(np.arange(-90,90,5.),labels=[True,False,False,True],fontsize=30)
proj.drawmeridians(np.arange(-180,180,5.),labels=[True,False,False,True],fontsize=30)
pc = proj.contourf(x2, y2, lev, range(75), cmap="jet")#,vmin=0., vmax=75.,zorder=3)
cbar = proj.colorbar(pc)
cbar.ax.tick_params(labelsize=30)
cbar.set_label('Level where stiffness is max [$m$]',size=40)
out_name ='lev_max_stiff.png'
plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1)
plt.close()

