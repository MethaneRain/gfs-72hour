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
        self.im_save_path =f"/Users/chowdahead/Desktop/{self.now.year}/{self.now.month:02d}_{self.now.day:02d}/"
        print(self.im_save_path)
            
        # Check to see if the folder already exists, if not create it
        if not os.path.isdir(self.im_save_path):
            os.makedirs(self.im_save_path)
            
        # Uncomment if you want to automatically change to the map folder    
        #os.chdir(im_save_path)
        
        # get current date and 00Z 
        self.start = datetime(self.now.year,self.now.month,self.now.day,0)
        
        # define time range you want the data for
        print(self.start)
        self.delta_t = 72
        self.end = self.start + timedelta(hours=self.delta_t)
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
        self.precip_tot_name = 'Total_precipitation_surface_Mixed_intervals_Accumulation'
        
        #self.query_list = [self.hgt_name,self.u_name,self.v_name]
        
        self.query_list = [self.mslp_name,self.u_src_name,self.v_src_name,self.precip_tot_name,
                          self.hgt_name,self.u_name,self.v_name]
        
        # Set the title font 
        self.title_font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 14,
        }
    def datetime_difference(self,date1,date2):
        '''
        Find time differences between two datetime instances
        ----------------------------------------------------
        
        Useful for finding FXXX number for file name and for title info for
        forecasted hour
        
        '''    
        diff = date2 - date1
        
        days, seconds = diff.days, diff.seconds
        hours = days * 24 + seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        
        return hours,minutes,seconds
    
    def get_data(self):
        # Request the GFS data from the thredds server
        gfs_url = f"https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/GFS_Global_0p25deg_{self.now.year}{self.now.month:02d}{self.now.day:02d}_0000.grib2/catalog.xml"
        gfs_cat = TDSCatalog(gfs_url)
        #https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'    
        dataset = list(gfs_cat.datasets.values())[0]
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
        print("done grabbing data!!\n-_-_-_-_-_-_-_-_-_-_-_\n")        
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
        
    
    
    def get_data_times(self,data,var):
        #var=precip_total_name
        '''
        Get all the forecast times in the data as strings
        -------------------------------------------------
        
        Arguments
        ---------
        
        
        Returns
        -------
        time_var: netCDF variable of all the times in the given ncss variable in dataset
        
        time_datetimes: list of datetime objects from ncss variable 
        
        time_strings: list of converted datetime values to strings
        * format %m/%d %H:%M
        '''
        # Get time into a datetime object
        time_var = data.variables[self.find_time_var(data.variables[var])]
        time_datetimes = num2date(time_var[:], time_var.units).tolist()
        print(type(time_datetimes),time_datetimes,"\n")
        time_strings = [t.strftime('%m/%d %H:%M') for t in time_datetimes]
            
        time_final = num2date(time_var[:].squeeze(), time_var.units)
        return time_var,time_strings,time_final
        
    def make_time_string(self,data,time_index,time_strings,time_final):
        '''
        Grab string of date and time for one specific time step in data.
        ie time_index=5 gets the sixth time of data. 
        
        Arguments
        ---------
        data: full dataset
        time_index: desired time step in dataset
        time_strings: list of requested dates in human readable string format
        time_final: list of requested datetimes
        
        
        Returns
        -------
        All returns are based off the requested single date index
        
        time_index: supplied time index
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
        time = time_strings[time_index]
        #print(f"string time: {time}")
        
        print(time_strings[time_index])
        
        # Set string for saved image file name
        file_time = str(time_final[0]).replace("-","_").replace(" ","_").replace(":","")[:-2]+"Z"
    
        # Set forecast date and hour  
        forecast_date = "{}".format(self.now.year)+'-'+str(time_strings[time_index]).replace("/","-")[:-5]
        forecast_hour = str(time_strings[time_index])[-5:]+"Z"
    
        # Set initialization date and hour 
        init_date = "{}".format(self.now.year)+'-'+str(time_strings[0]).replace("/","-")[:-5]
        init_hour = str(time_strings[0]).replace("/","-")[-5:]+"Z"
        
        print(f"time index: {time_index},\n\
            string time: {time},\n\
            filename time: {file_time},\n\
            forecast date: {forecast_date},\n\
            forecast hour: {forecast_hour},\n\
            initialization date: {init_date},\n\
            initialization hour: {init_hour}\n")
        return time,file_time,forecast_date,forecast_hour,init_date,init_hour
    
    
    def get_var_isobaric_num(self,data,var):
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
            print("no isobaricX variable name!")
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




    
    
    
    
    
    
    
    
    
    
    
   