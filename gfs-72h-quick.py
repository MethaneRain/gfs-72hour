#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 18:44:37 2020

@author: Justin Richling
"""


# System Tools
import os,sys
            
# Importing Datetime Libraries
from datetime import datetime, timedelta
        
# MetPy
from metpy.units import units
            
# CartoPy Map Plotting Libraires
import cartopy.crs as ccrs
import cartopy.feature as cfeature
            
# Numerical and Scientific Libraries
import numpy as np
from scipy.ndimage import gaussian_filter
            
# Accessing Data from XLM Catalog via Siphon Libraries
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
            
# MetPy Libraries
from metpy.plots import add_metpy_logo
            
# NetCDF Libraries
from netCDF4 import num2date
            
# Matplotlib Plotting Libraries
import matplotlib.pyplot as plt
from matplotlib import patheffects
            
# Warnings
import warnings
warnings.filterwarnings('ignore')
            

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
class GFS_72hour_Maps:
    '''
    Class to collect all my 72-hour GFS forecast maps. 
    
    Data comes from THREDDS via Siphon package and MetPy, much thanks
    to the devs and everyoone involved in these packages and forums.
    
    Variables Plotted:
    -----------------
    
    * Potential Vorticity by Pressure and 500mb Heights
    * 250mb Jets
    * MSLP, Hi/Lows, and 1000-500mb Thickness
    * 500mb Heights and Absolute Vorticity
    * 
    '''
    def __init__(self):
        
        self.plotcrs = ccrs.PlateCarree()    


        self.now = datetime.utcnow()
        self.im_save_path =f"/Users/chowdahead/Desktop/Weather_Blog/{self.now.year}/{self.now.month}_{self.now.day}/"
        print(self.im_save_path)
            
        # Check to see if the folder already exists, if not create it
        if not os.path.isdir(self.im_save_path):
            os.makedirs(self.im_save_path)
            
        # Uncomment if you want to automatically change to the map folder    
        #os.chdir(im_save_path)
        
        # get current date and time
        #now = forecast_times[0]
        self.start = datetime(self.now.year,self.now.month,self.now.day,0)
        # define time range you want the data for
        print(self.start)
        delta_t = 72
        self.end = self.start + timedelta(hours=delta_t)
        self.extent = [-130,-60,20,60]
                
        return
    
    def get_data(self):
        # Request the GFS data from the thredds server
        gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml')
            
        dataset = list(gfs.datasets.values())[1]
        #print(dataset.access_urls)
            
        # Create NCSS object to access the NetcdfSubset
        ncss = NCSS(dataset.access_urls['NetcdfSubset'])
        
        # query the data from the server
        query = ncss.query()
        query.time_range(self.start, self.end)
        query.lonlat_box(north=80, south=0, east=310, west=200)
        query.accept('netcdf4')
        data = ncss.get_data(query)
                
        return data
    
    
    def find_time_var(var, time_basename='time'):
        '''
        Thanks to the crew over at Metpy for this handy little function
        
        Grabs time listed in variable's metadata
        '''
        for coord_name in var.coordinates.split():
            if coord_name.startswith(time_basename):
                return coord_name
        raise ValueError('No time variable found for ' + var.name)
        
    
    
        
     
    def get_time_string(self,data,time_index):
        # File and Title Times
        #---------------------------------------------------------------------------------------------------
        
        # Get time into a datetime object
        time_var = data.variables[self.find_time_var(data.variables["MSLP_Eta_model_reduction_msl"])]
        time_var = num2date(time_var[:], time_var.units).tolist()
        time_strings = [t.strftime('%m/%d %H:%M') for t in time_var]
            
        time_var = data.variables[self.find_time_var(data.variables["MSLP_Eta_model_reduction_msl"])]
        time_final = num2date(time_var[:].squeeze(), time_var.units)
        
        # Time index for data variables
        time = time_strings[time_index]
        print(time)
    
        # Set string for saved image file name
        file_time = str(time_final[0]).replace("-","_").replace(" ","_").replace(":","")[:-2]+"Z"
    
        # Set forecast date and hour  
        forecast_date = "{}".format(self.now.year)+'-'+time_strings[time_index].replace("/","-")[:-5]
        forecast_hour = time_strings[time_index][-5:]+"Z"
    
        # Set initialization date and hour 
        init_date = "{}".format(self.now.year)+'-'+time_strings[0].replace("/","-")[:-5]
        init_hour = time_strings[0].replace("/","-")[-5:]+"Z"
        
        return time,file_time,forecast_date,forecast_hour,init_date,init_hour,time_strings
    
    
    
    
    def get_var_isobaric_num(data,var):
        '''
        thredds use different isobaric variable names depending on the variable,
        ie, isobaric1, isobaric7, etc. Thus grab the isobaric number variable name
        
        input variable that has an isobaricX level, if the variable does have an isobaric number it will be set
        '''
        iso = data.variables[var].coordinates[:]
        finder = iso.find("isobaric")
        if "isobaric" in iso:
            #print(iso[finder:finder+9])
            iso_num = iso[finder:finder+9]
        else:
            print("no isobaricX varible name!")
        type(iso_num),iso_num
        return iso_num   
    

    def make_map(self):
        
        # Setup Figure
        #---------------------------------------------------------------------------------------------------    
        fig = plt.figure(figsize=(17., 11.))
            
        #add_metpy_logo(fig, 25, 950, size='small')
            
        # Add the Map 
        #---------------------------------------------------------------------------------------------------
        ax = plt.subplot(111, projection=self.plotcrs)
            
        # Set extent and plot map lines
        ax.set_extent(self.extent)
    
        ax.coastlines(resolution='50m')
            
        # Add State/Country Boundaries to Plot
        #---------------------------------------------------------------------------------------------------    
        state_borders = cfeature.NaturalEarthFeature(
            category='cultural', name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
        ax.add_feature(state_borders, edgecolor='k', linewidth=1, zorder=3,linestyle="--")
    
        country_borders = cfeature.NaturalEarthFeature(category='cultural',
            name='admin_0_countries',scale='50m', facecolor='none')
        ax.add_feature(country_borders,edgecolor='k',linewidth=0.9)
    
        
        return fig,ax




    
    
    
    
    
    
    
    
    
    
    
   