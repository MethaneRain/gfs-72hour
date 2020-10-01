#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 18:00:30 2020

@author: Justin Richling
"""

# System Tools
import os

# MetPy
from metpy.units import units

# Numerical and Scientific Libraries
import numpy as np
from scipy.ndimage import gaussian_filter

from gfs-72h-quick import GFS_72hour_Maps
data = GFS_72hour_Maps.get_data()

make_map = GFS_72hour_Maps.make_map()


vort_name = "Absolute_vorticity_isobaric"
hgt_name = "Geopotential_height_isobaric"
u_src_name = "u-component_of_wind_height_above_ground"
v_src_name = "v-component_of_wind_height_above_ground"
sfc_gust_name = 'Wind_speed_gust_surface'
pv_press_name = "Pressure_potential_vorticity_surface"
mslp_name = "MSLP_Eta_model_reduction_msl"
upflux_rad_name = "Upward_Long-Wave_Radp_Flux_atmosphere_top_Mixed_intervals_Average"
u_name = 'u-component_of_wind_isobaric'
v_name = 'v-component_of_wind_isobaric'


vort = data.variables[vort_name][:]
hgt = data.variables[hgt_name][:]

pv = data.variables[pv_press_name][:]
mslp = data.variables[mslp_name][:]
upflux_rad = data.variables[upflux_rad_name][:]

u = data.variables[u_name][0] * units('m/s')
v = data.variables[v_name][0] * units('m/s')
lev_250 = np.where(data.variables['isobaric'][:] == 25000)[0][0]
u_250 = data.variables[u_name][:, lev_250, :, :]
v_250 = data.variables[v_name][:, lev_250, :, :]

u_sfc = data.variables[u_src_name][:] * units('m/s')
v_sfc = data.variables[v_src_name][:] * units('m/s')


# Grab pressure levels
plev = list(data.variables['isobaric'][:])
# Grab pressure level data
hght_1000 = data.variables['Geopotential_height_isobaric'][:, plev.index(1000)]
hght_500 = data.variables['Geopotential_height_isobaric'][:, plev.index(500)]
# Calculate and smooth 1000-500 hPa thickness
thickness_1000_500 = gaussian_filter(hght_500 - hght_1000, sigma=3.0)

# Pull out the lat and lon data
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
    
# Combine 1D latitude and longitudes into a 2D grid of locations
# use this for the High/Low function 
lon_2d, lat_2d = np.meshgrid(lons, lats)

            
os.chdir("map-scripts/")
import HiLo_map as hilo
hilo.HiLo_thickness_map(time_index,extent,im_save_path,get_time_string,title_font,make_map,
                       lons,lats,mslp)
