#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 18:00:30 2020

@author: Justin Richling
"""

# System Tools
import os

from gfs_72h_quick import GFS_72hour_Maps

gfs = GFS_72hour_Maps()

data = gfs.get_data()

fig,ax = gfs.make_map()

#gfs.g

'''
vort = data.variables[gfs.vort_name][:]
hgt = data.variables[gfs.hgt_name][:]

pv = data.variables[gfs.pv_press_name][:]
upflux_rad = data.variables[gfs.upflux_rad_name][:]

u = data.variables[gfs.u_name][0] * units('m/s')
v = data.variables[gfs.v_name][0] * units('m/s')
lev_250 = np.where(data.variables['isobaric'][:] == 25000)[0][0]
u_250 = data.variables[gfs.u_name][:, lev_250, :, :]
v_250 = data.variables[gfs.v_name][:, lev_250, :, :]

u_sfc = data.variables[gfs.u_src_name][:] * units('m/s')
v_sfc = data.variables[gfs.v_src_name][:] * units('m/s')


# Grab pressure levels
plev = list(data.variables['isobaric'][:])
# Grab pressure level data
hght_1000 = data.variables['Geopotential_height_isobaric'][:, plev.index(1000)]
hght_500 = data.variables['Geopotential_height_isobaric'][:, plev.index(500)]
# Calculate and smooth 1000-500 hPa thickness
thickness_1000_500 = gaussian_filter(hght_500 - hght_1000, sigma=3.0)
'''

#mslp = data.variables[gfs.mslp_name][:]
    
# Combine 1D latitude and longitudes into a 2D grid of locations
# use this for the High/Low function 
#lon_2d, lat_2d = np.meshgrid(lons, lats)

#extent = gfs.extent     
#im_save_path = gfs.im_save_path
#im_save_path = GFS_72hour_Maps.
#get_time_string = gfs.get_time_string(data,5)

print("\n\nchanging dirs for HiLo script...\n")     
os.chdir("/Users/chowdahead/Documents/gfs-72hour/map-scripts")
print(os.getcwd())
print(os.listdir())
print("importing HiLo script...\n") 
import sys
sys.path.append('/Users/chowdahead/Documents/gfs-72hour/map-scripts')
import HiLo_map as hilo
print("import done...\ntrying map plot function...\n")
hilo.HiLo_thickness_map(data,
                       5,
                       fig,
                       ax)
fig
