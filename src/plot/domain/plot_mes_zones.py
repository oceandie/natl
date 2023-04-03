#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

config = 'orca025'

if config == 'orca025':
   #MEs = "/data/users/dbruciaf/GS/orca025/MEs_GO8/polyg_1/4env_MEs_000_015_020/bathymetry.MEs_4env_000-015-020_v1_maxdep_3750.0.nc"
   #MEs = "/data/users/dbruciaf/GS/orca025/MEs_GO8/polyg_2/bathymetry.loc_area.dep4500_pol2_sig4_itr1.nc" 
    MEs = "/net/home/h01/dbruciaf/mod_dev/local-ME_GS/src/loc_area/bathymetry.loc_area.dep4350_pol2_sig4_itr1.nc"
    MEs = "/home/h01/dbruciaf/mod_dev/local-ME_GS/src/loc_area/bathymetry.loc_area.dep0.002_pol3_sig2_itr1.nc"
else:
   MEs = "/home/h01/dbruciaf/mod_dev/local-ME_GS/src/envelopes/orca12/bathymetry.MEs_4env_020-015-015_maxdep_3750.0.nc"

nc       = nc4.Dataset(MEs,"r")
bathy   = np.array(nc.variables["Bathymetry"])[:,:]
lon   = np.array(nc.variables["nav_lon"])[:,:]
lat   = np.array(nc.variables["nav_lat"])[:,:]
mes_gs = np.array(nc.variables["loc_area"])[:,:]
s2z_gs = np.array(nc.variables["s2z_msk"])[:,:]
nc.close()
land = np.ma.masked_where(bathy>0, bathy*0.+1. )
bathy[s2z_gs==0]=np.nan
mes_gs = np.ma.array(mes_gs)
mes_gs = np.ma.masked_where(mes_gs==0, mes_gs)
s2z_gs = np.ma.array(s2z_gs)
s2z_gs = np.ma.masked_where(np.logical_or(s2z_gs==0,s2z_gs==2), s2z_gs)

LLcrnrlon = -90. 
LLcrnrlat =  20. 
URcrnrlon = -30.  
URcrnrlat =  70.

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


cn_lev = [250., 500., 1000.]
proj.contour(x2, y2, bathy, levels=cn_lev,colors='black',linewidths=1.5,zorder=4)
proj.contour(x2, y2, bathy, levels=[4350.],colors='black',linewidths=4.,zorder=4)

proj.pcolor(x2, y2, mes_gs, cmap = 'autumn',zorder=3,alpha=0.5)
proj.pcolor(x2, y2, s2z_gs, cmap = 'winter_r',zorder=2,alpha=0.5)
#proj.pcolor(x2, y2, land, cmap = 'binary_r',zorder=4)

#proj.scatter(x2[960,958] , y2[960,958], s=100, color='black')

out_name =config+'_mes_areas_natl.png'
plt.savefig(out_name,bbox_inches="tight", pad_inches=0.1)
plt.close()

