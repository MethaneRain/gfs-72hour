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
        
        self.plotcrs = ccrs.LambertConformal()    


        self.now = datetime.utcnow()
        self.im_save_path =f"{self.now.year}/{self.now.month}_{self.now.day}/"
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
        self.vort_name = "Absolute_vorticity_isobaric"
        self.hgt_name = "Geopotential_height_isobaric"
        self.u_src_name = "u-component_of_wind_height_above_ground"
        self.v_src_name = "v-component_of_wind_height_above_ground"
        self.sfc_gust_name = 'Wind_speed_gust_surface'
        self.pv_press_name = "Pressure_potential_vorticity_surface"
        self.mslp_name = "MSLP_Eta_model_reduction_msl"
        self.upflux_rad_name = "Upward_Long-Wave_Radp_Flux_atmosphere_top_Mixed_intervals_Average"
        self.u_name = 'u-component_of_wind_isobaric'
        self.v_name = 'v-component_of_wind_isobaric'
        
        self.query_list = [self.mslp_name,self.u_src_name,self.v_src_name]
        
        # Set the title font 
        self.title_font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 14,
        }
    
    def get_data(self):
        # Request the GFS data from the thredds server
        gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml')
            
        dataset = list(gfs.datasets.values())[1]
        #print(dataset.access_urls)
            
        # Create NCSS object to access the NetcdfSubset
        ncss = NCSS(dataset.access_urls['NetcdfSubset'])
        
        # query the data from the server
        query = ncss.query()
        query.time_range(self.start,self.end)
        query.lonlat_box(north=80, south=0, east=310, west=200)
        query.accept('netcdf4')
        
        print("-----------------------------------------\n"\
              +"Sit back....\nOr get your coffee....\nOr do a Sudoku....\n"\
              +"-----------------------------------------\n")
        print("qeueing data...") 
        #query.variables(str(self.query_list)).add_lonlat(True)
        for i in self.query_list:
            query.variables(i) 
        #query.variables(vort_name,hgt_name,pv_press_name,mslp_name,upflux_rad_name,u_name,v_name,
        #               u_src_name,v_src_name,sfc_gust_name).add_lonlat(True)
        print("\ndone qeueing data.\n\ngrabbing data...\n")  
            
        # Request data for the variables you want to use
        self.data = ncss.get_data(query)
        print("done grabbing data!!\n-_-_-_-_-_-_-_-_-_-_-_")        
        return self.data
    
    
    def find_time_var(self,var, time_basename='time'):
        '''
        Thanks to the crew over at Metpy for this handy little function
        
        Grabs time listed in variable's metadata
        '''
        for coord_name in var.coordinates.split():
            if coord_name.startswith(time_basename):
                return coord_name
        raise ValueError('No time variable found for ' + var.name)
        
    
    
    def get_data_times(self,data):
        '''
        Get all the forecast times in the data as strings
        -------------------------------------------------
        '''
        # Get time into a datetime object
        time_var = data.variables[self.find_time_var(data.variables["MSLP_Eta_model_reduction_msl"])]
        time_var = num2date(time_var[:], time_var.units).tolist()
        self.time_strings = [t.strftime('%m/%d %H:%M') for t in time_var]
            
        time_var = data.variables[self.find_time_var(data.variables["MSLP_Eta_model_reduction_msl"])]
        self.time_final = num2date(time_var[:].squeeze(), time_var.units)
        
        return self.time_strings,self.time_final
        
    def get_time_string(self,data,time_index):
        '''
        Grab string of date and time for one specific time step in data.
        ie time_index=5 gets the sixth time of data. 
        
        Arguments
        ---------
        data: full dataset
        time_index: desired time step in dataset
        
        Returns
        -------
        time_index: time index
        time: string time
        file_time: filename time
        forecast_date: forecast date
        forecast_hour: forecast hour
        init_date: initialization date
        init_hour: initialization hour
        
        '''
        # File and Title Times
        #---------------------------------------------------------------------------------------------------
        
        # Time index for data variables
        time = self.time_strings[time_index]
        #print(f"string time: {time}")
    
        # Set string for saved image file name
        file_time = str(self.time_final[0]).replace("-","_").replace(" ","_").replace(":","")[:-2]+"Z"
    
        # Set forecast date and hour  
        forecast_date = "{}".format(self.now.year)+'-'+self.time_strings[time_index].replace("/","-")[:-5]
        forecast_hour = self.time_strings[time_index][-5:]+"Z"
    
        # Set initialization date and hour 
        init_date = "{}".format(self.now.year)+'-'+self.time_strings[0].replace("/","-")[:-5]
        init_hour = self.time_strings[0].replace("/","-")[-5:]+"Z"
        print(f"time index: {time_index},\n\
            string time: {time},\n\
            filename time: {file_time},\n\
            forecast date: {forecast_date},\n\
            forecast hour: {forecast_hour},\n\
            initialization date: {init_date},\n\
            initialization hour: {init_hour}\n")
        return time,file_time,forecast_date,forecast_hour,init_date,init_hour
    
    
    
    @staticmethod
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




    
    
    
    
    
    
    
    
    
    
    
   