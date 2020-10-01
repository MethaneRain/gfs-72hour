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
        '''
        Default values
        
        * date -> current date
        * initialization hour -> 0000Z
        * image save path -> cwd
        * countour levels -> 100
        * extent -> CONUS
        * contour color -> black
        '''
        
        
        
        
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
        self.map_dict = {'vort':{"ax_title":'GFS 0.5$^{o}$\n500mb Heights (m) and Abs Vorticity '+ "("+"$\mathregular{s^{-1}}$"+")",
                            "long_name":'Absolute_vorticity_isobaric',
                            "save_path":"GFS/Vorticity/",
                            "im_name":"GFS_0p5_Vort_Heights_500mb"},
                    'mslp':{"ax_title":'GFS 0.5$^{o}$\nMean Sea Level Pressure (hPa)',
                            "long_name":"MSLP_Eta_model_reduction_msl",
                            "save_path":"GFS/MSLP/",
                            "im_name":"GFS_0p5_mslp"},
                    'jets':{"ax_title":'GFS 0.5$^{o}$\n500hPa Heights (m) and 250hPa Windspeed (m/s)',
                            "long_name":"",
                            "save_path":"GFS/Wind/",
                            "im_name":"GFS_0p5_WindSpeed"},
                    'pv':{"ax_title":"",
                            "long_name":"",
                            "save_path":"",
                            "im_name":""},
                }
        
        
        
        # Set the title font 
        self.title_font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 14,
        }
        self.now = datetime.utcnow()
        self.now_start = datetime(self.now.year,self.now.month,self.now.day,0,0)
        # define time range you want the data for
        self.start = self.now_start
        print(self.start)
        self.delta_t = 72
        self.end = self.now_start + timedelta(hours=self.delta_t)
            
        # Set a path to save the plots with string format for the date to set the month and day
        self.im_save_path ="/path/to/saved/images/"
        self.im_save_path =f"/Users/chowdahead/Desktop/Weather_Blog/{self.now.year}/{self.now.month}_{self.now.day}/"
        print(self.im_save_path)
            
        # Check to see if the folder already exists, if not create it
        if not os.path.isdir(self.im_save_path):
            os.makedirs(self.im_save_path)
            
        # Uncomment if you want to automatically change to the map folder    
        #os.chdir(im_save_path)
        
        # get current date and time
        #now = forecast_times[0]
        
        # initialization hour for model forecast
        self.init_hour = "0000"
        
        # current working drive for saved images
        self.im_save_path = "./"
        
        # 
        self.arg = ""
        
        # contour levels for plotting
        self.clevs = 100
        
        # map lat/lon extent
        self.extent = [-130., -70, 20., 60.]
        # Set Projection of Data
        self.datacrs = ccrs.PlateCarree()
            
        # Set Projection of Plot
        self.plotcrs = ccrs.LambertConformal(central_latitude=[30, 60], central_longitude=-100)
            
        # Add Map Features
        self.states_provinces = cfeature.NaturalEarthFeature(category='cultural',
        name='admin_1_states_provinces_lakes',scale='50m', facecolor='none')
            
        self.country_borders = cfeature.NaturalEarthFeature(category='cultural',
        name='admin_0_countries',scale='50m', facecolor='none')
            
        # Colorbar Axis Placement (under figure)
        self.colorbar_axis = [0.183, 0.09, 0.659, 0.03] # [left, bottom, width, height]
            
         
        
        # plot contour line color
        self.colors = "k"
        
        # Request the GFS data from the thredds server
        gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml')
            
        self.dataset = list(gfs.datasets.values())[1]
        print(self.dataset.access_urls)
            
        # Create NCSS object to access the NetcdfSubset
        self.ncss = NCSS(self.dataset.access_urls['NetcdfSubset'])
        
        # query the data from the server
        self.query = self.ncss.query()
        self.query.time_range(self.start, self.end)
        self.north = 80
        self.south = 0
        self.east = 310
        self.west = 200
        self.Lat = self.query.lonlat_box(north=self.north, south=self.south, east=self.east, west=self.west)
        
        self.query.variables("MSLP_Eta_model_reduction_msl").add_lonlat(True)
        # Request data for the variables you want to use
        self.data_o = self.ncss.get_data(self.query)
        #self.get_time_string(self,self.now.year,time_index)
        
        # Pull out the lat and lon data
        self.lats = self.data_o.variables['lat'][:]
        self.lons = self.data_o.variables['lon'][:]
            
        # Combine 1D latitude and longitudes into a 2D grid of locations
        # use this for the High/Low function 
        self.lon_2d, self.lat_2d = np.meshgrid(self.lons, self.lats)
                
        self.mslp = self.data_o.variables[self.mslp_name][:]
                
    
    def find_time_var(self,var, time_basename='time'):
        '''
        Thanks to the crew over at Metpy for this handy little function
        
        Grabs time listed in variable's metadata
        '''
        for coord_name in var.coordinates.split():
            if coord_name.startswith(time_basename):
                return coord_name
        raise ValueError('No time variable found for ' + var.name)
        
    
    def get_variables(self,var_name_list):
        '''
        get and/or change dataset from thredds
        
        -------
        Args:
        * <em>list</em> of strings corresponding to variable names in THREDDS
        
        '''
       
        self.query = self.ncss.query()
        for var in var_name_list:
            self.query.variables(var).add_lonlat(True)
        
        print("\ndone queing data...\n\ngrabbing data...\n")  
            
        # Request data for the variables you want to use
        self.data = self.ncss.get_data(self.query)
        print(type(self.data))
        print("done grabbing data...\n\ngrabbing variables and times...\n")
        return self.data
        # Get time into a datetime object
        
     
    def get_time_string(self,time_index):
        # File and Title Times
        #---------------------------------------------------------------------------------------------------
        
        # Get time into a datetime object
        time_var = self.data_o.variables[self.find_time_var(self.data_o.variables["MSLP_Eta_model_reduction_msl"])]
        time_var = num2date(time_var[:], time_var.units).tolist()
        time_strings = [t.strftime('%m/%d %H:%M') for t in time_var]
            
        time_var = self.data_o.variables[self.find_time_var(self.data_o.variables["MSLP_Eta_model_reduction_msl"])]
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
    
    
    
    
    def get_var_isobaric_num(self,var):
        '''
        thredds use different isobaric variable names depending on the variable,
        ie, isobaric1, isobaric7, etc. Thus grab the isobaric number variable name
        
        input variable that has an isobaricX level, if the variable does have an isobaric number it will be set
        '''
        iso = self.data.variables[var].coordinates[:]
        finder = iso.find("isobaric")
        if "isobaric" in iso:
            #print(iso[finder:finder+9])
            iso_num = iso[finder:finder+9]
        else:
            print("no isobaricX varible name!")
        type(iso_num),iso_num
        return iso_num   
    

    def make_map(self, extent):
        
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
            
    
    def HiLo_thickness_map(self,time_index,thickness_plot,u_sfc,v_sfc):
        '''
        Method to plot the 500mb heights and absolute vorticity
        -------------------------------------------------------
    
        argument: time index - for GFS forecast hour interval
        '''
        thickness_plot=thickness_plot
        
        
        # Setup Contour Label Options
        #---------------------------------------------------------------------------------------------------    
        kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
        'rightside_up': True, 'use_clabeltext': True}
        
        fig,ax = self.make_map(self.extent)
        
        # Plot Title
        #---------------------------------------------------------------------------------------------------
        
        time,file_time,forecast_date,forecast_hour,init_date,init_hour,time_strings = self.get_time_string(time_index)
        ax.set_title('GFS 0.5$^{o}$\n500hPa Heights (m) and PVU (hPa)', 
        size=10, loc='left',fontdict=self.title_font)
    
        ax.set_title(f"Init Hour: {init_date} {init_hour}\nForecast Hour: {forecast_date} {forecast_hour}",
        size=10, loc='right',fontdict=self.title_font)
            
        ax.stock_img()
        # 250hPa Jet
        #---------------------------------------------------------------------------------------------------
        
        #u_sfc = self.data.variables[self.u_src_name][self.time_strings.index(time),0] * units('m/s')
        #print(u_sfc.shape)
        #v_sfc = self.data.variables[self.v_src_name][self.time_strings.index(time),0] * units('m/s')
        
        
        # Grab pressure levels
        #plev = list(self.data.variables['isobaric'][:])
            
        # Plot MSLP
        clevmslp = np.arange(800., 1120., 4)
        mslp_smooth = gaussian_filter(self.mslp[time_strings.index(time),:,:],sigma=3.0)
        cs2 = ax.contour(self.lons, self.lats, mslp_smooth/100., clevmslp, colors='k', linewidths=1.25,
                         linestyles='solid', transform=ccrs.PlateCarree())
        clbls = plt.clabel(cs2, **kw_clabels)
        plt.setp(clbls, path_effects=[
            patheffects.withStroke(linewidth=3, foreground="w")])
    
        # Use definition to plot H/L symbols
        #plot_maxmin_points(ax,self.lon_2d, self.lat_2d, mslp_smooth/100., 'max', 50, symbol='H', color='b',  transform=ccrs.PlateCarree())
        #plot_maxmin_points(ax,self.lon_2d, self.lat_2d, mslp_smooth/100., 'min', 25, symbol='L', color='r', transform=ccrs.PlateCarree())
        
        
        lon_slice = slice(None, None, 15)
        lat_slice = slice(None, None, 15)
        ax.barbs(self.lons[lon_slice], self.lats[lat_slice],
                 u_sfc[lon_slice, lat_slice],
                 v_sfc[lon_slice, lat_slice],
                 np.sqrt(u_sfc[lon_slice, lat_slice]**2+v_sfc[lon_slice, lat_slice]**2),
                 cmap="magma",
                 #isen_u[lon_slice, lat_slice].to('knots').magnitude,
             #isen_v[lon_slice, lat_slice].to('knots').magnitude,
             transform=ccrs.PlateCarree(), zorder=20,length=6,
             regrid_shape=20,)
            #sizes=dict(spacing=0.2,width=0.3)) # height=8
        
        '''
        ax.barbs(lons[lon_slice], lats[lat_slice],
                 u_sfc[lon_slice, lat_slice],
                 v_sfc[lon_slice, lat_slice],
                 np.sqrt(u_sfc[lon_slice, lat_slice]**2+v_sfc[lon_slice, lat_slice]**2)
                 cmap="jet",
                 #isen_u[lon_slice, lat_slice].to('knots').magnitude,
             #isen_v[lon_slice, lat_slice].to('knots').magnitude,
             transform=ccrs.PlateCarree(), zorder=20,length=6.75,
             regrid_shape=20,barbcolor="w",color="w")
        '''
        # Save Figure
        #---------------------------------------------------------------------------------------------------    
        HILO = self.im_save_path+"GFS/HiLo/"
        if not os.path.isdir(HILO):
            os.makedirs(HILO)
            
        time_index *= 3
        if time_index < 10:
            times = f"0{time_index}"
        else:
           times = f"{time_index}"
        
        
        #if thickness_plot == True:
        '''
            # Grab pressure level data
            hght_1000 = self.data.variables['Geopotential_height_isobaric'][self.time_strings.index(time), plev.index(1000)]
            hght_500 = self.data.variables['Geopotential_height_isobaric'][self.time_strings.index(time), plev.index(500)]
            # Calculate and smooth 1000-500 hPa thickness
            thickness_1000_500 = gaussian_filter(hght_500 - hght_1000, sigma=3.0)
            # Plot thickness with multiple colors
            clevs = (np.arange(0, 5400, 60),
                     np.array([5400]),
                     np.arange(5460, 7000, 60))
            colors = ('tab:blue', 'b', 'tab:red')
            kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                          'rightside_up': True, 'use_clabeltext': True}
            for clevthick, color in zip(clevs, colors):
                cs = ax.contour(self.lons, self.lats, thickness_1000_500, levels=clevthick, colors=color,
                                linewidths=1.0, linestyles='dashed', transform=ccrs.PlateCarree())
                plt.clabel(cs, **kw_clabels)
        '''
       #     outfile = f"{HILO}GFS_0p5_HiLo_thickness_{file_time}_F{times}.png"
       #else:
       #     outfile = f"{HILO}GFS_0p5_HiLo_{file_time}_F{times}.png"
        
        
        outfile = f"{HILO}GFS_0p5_HiLo_thickness_{file_time}_F{times}.png"
        fig.savefig(outfile,bbox_inches='tight',dpi=120)
        #plt.close(fig)
        #print(time,file_time,forecast_date,forecast_hour,init_date,init_hour)
        
        #return file_time
    
    print("ahhhh")
    
    
    
    
    
    
    
    
    
    
    
   