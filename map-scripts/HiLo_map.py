#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 00:30:10 2020

@author: chowdahead
"""

def HiLo_thickness_map(time_index,extent,im_save_path,get_time_string,title_font,make_map,
                       lons,lats,mslp):
    '''
    Method to plot the 500mb heights and absolute vorticity
    -------------------------------------------------------
    
    argument: time index - for GFS forecast hour interval
    '''
    thickness_plot=False
    import os
    
    # MetPy
    from metpy.units import units
    
    # Numerical and Scientific Libraries
    import numpy as np
    from scipy.ndimage import gaussian_filter
    
    import matplotlib.pyplot as plt
    
    # CartoPy Map Plotting Libraires
    import cartopy.crs as ccrs
    
    #from HiLo import plot_maxmin_points as plot_maxmin_points
    
    from matplotlib import patheffects
    
    
    # Setup Contour Label Options
    #---------------------------------------------------------------------------------------------------    
    kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
    'rightside_up': True, 'use_clabeltext': True}
    
    fig,ax = make_map()
    
    # Plot Title
    #---------------------------------------------------------------------------------------------------
    
    time,file_time,forecast_date,forecast_hour,init_date,init_hour,time_strings = get_time_string()
    ax.set_title('GFS 0.5$^{o}$\n500hPa Heights (m) and PVU (hPa)', 
    size=10, loc='left',fontdict=title_font)
    
    ax.set_title(f"Init Hour: {init_date} {init_hour}\nForecast Hour: {forecast_date} {forecast_hour}",
    size=10, loc='right',fontdict=title_font)
        
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
    mslp_smooth = gaussian_filter(mslp[time_strings.index(time),:,:],sigma=3.0)
    cs2 = ax.contour(lons, lats, mslp_smooth/100., clevmslp, colors='k', linewidths=1.25,
     linestyles='solid', transform=ccrs.PlateCarree())
    clbls = plt.clabel(cs2, **kw_clabels)
    plt.setp(clbls, path_effects=[
        patheffects.withStroke(linewidth=3, foreground="w")])
    
    # Use definition to plot H/L symbols
    #plot_maxmin_points(ax,self.lon_2d, self.lat_2d, mslp_smooth/100., 'max', 50, symbol='H', color='b',  transform=ccrs.PlateCarree())
    #plot_maxmin_points(ax,self.lon_2d, self.lat_2d, mslp_smooth/100., 'min', 25, symbol='L', color='r', transform=ccrs.PlateCarree())
    
    '''
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
    HILO = im_save_path+"GFS/HiLo/"
    if not os.path.isdir(HILO):
        os.makedirs(HILO)
        
    time_index *= 3
    if time_index < 10:
        times = f"0{time_index}"
    else:
        times = f"{time_index}"
    
    
    if thickness_plot == True:
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
        outfile = f"{HILO}GFS_0p5_HiLo_thickness_{file_time}_F{times}.png"
    else:
        outfile = f"{HILO}GFS_0p5_HiLo_{file_time}_F{times}.png"
    
    
    #outfile = f"{HILO}GFS_0p5_HiLo_thickness_{file_time}_F{times}.png"
    fig.savefig(outfile,bbox_inches='tight',dpi=120)
    plt.close(fig)
    print(time,file_time,forecast_date,forecast_hour,init_date,init_hour)
    
    return file_time