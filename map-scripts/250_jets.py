#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 18:18:52 2020

@author: Justin Richling
"""
    
def jet250_map(data,
                       time_index,
                       fig,
                       ax):
    '''
    Method to plot the 250mb jets and 500mb heights
    -------------------------------------------------------
    
    argument: time index - for GFS forecast hour interval
    '''
   
    import os
        
    # Numerical and Scientific Libraries
    import numpy as np
    from scipy.ndimage import gaussian_filter
    
    import matplotlib.pyplot as plt
    
    # CartoPy Map Plotting Libraires
    import cartopy.crs as ccrs
    
    #from HiLo import plot_maxmin_points as plot_maxmin_points
    
    from matplotlib import patheffects
    
    # MetPy
    from metpy.units import units
    
     
    os.chdir("../")
    from gfs_72h_quick import GFS_72hour_Maps
    
    gfs = GFS_72hour_Maps()
    time_strings,_ = gfs.get_data_times(data)
    title_font = gfs.title_font
    ref_cmap = "magma"
    
    u_name = gfs.u_name
    v_name = gfs.v_name
    hgt_name = gfs.hgt_name
    hgt = data.variables[hgt_name][:]
    #u = data.variables[u_name][0] * units('m/s')
    #v = data.variables[v_name][0] * units('m/s')
    lev_250 = np.where(data.variables['isobaric'][:] == 25000)[0][0]
    u_250 = data.variables[u_name][:, lev_250, :, :]
    v_250 = data.variables[v_name][:, lev_250, :, :]

    
    
    #mslp_name = gfs.mslp_name
    im_save_path = gfs.im_save_path
    
    print(list(data.variables))
    lons = data.variables["lon"][:]
    lats = data.variables["lat"][:]    
    # Setup Contour Label Options
    #---------------------------------------------------------------------------------------------------    
    kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
    'rightside_up': True, 'use_clabeltext': True}
        
    # Setup Figure
    #---------------------------------------------------------------------------------------------------    
    # Plot Title
    #---------------------------------------------------------------------------------------------------
    time,file_time,forecast_date,forecast_hour,init_date,init_hour = gfs.get_time_string(data,time_index)
    #time,file_time,forecast_date,forecast_hour,init_date,init_hour,time_strings = get_time_string()


    # Plot Title
    #---------------------------------------------------------------------------------------------------
    ax.set_title('GFS 0.5$^{o}$\n500hPa Heights (m) and 250hPa Windspeed (m/s)', 
    size=10, loc='left',fontdict=title_font)

    ax.set_title(f"Init Hour: {init_date} {init_hour}\nForecast Hour: {forecast_date} {forecast_hour}",
    size=10, loc='right',fontdict=title_font)
        
    
    # 300hPa Jet
    #---------------------------------------------------------------------------------------------------
    
    jet_levs = np.arange(np.sqrt(u_250**2+v_250**2).min(),np.sqrt(u_250**2+v_250**2).max(),2)
    
    jet = np.sqrt(u_250**2+v_250**2)
    jet2 = np.ma.masked_where(jet<20,jet)
    jet2_2 = np.ma.masked_where(jet>40,jet)
    
    
    cs2 = ax.contourf(lons, lats, jet2,jet_levs,vmin=20,
    transform=ccrs.PlateCarree(),cmap=ref_cmap)
    
    #cs2_2 = ax.contourf(lons, lats, jet2_2,jet_levs,
    #transform=ccrs.PlateCarree(),cmap='magma',alpha=0.2)
    
    cbar = plt.colorbar(cs2,orientation="horizontal") #,ticks=ticks ticks=range(998,1030,8)
    posn = ax.get_position()
    print(posn.x0,posn.x1)
    cbar.ax.set_position([posn.x0+0.001, posn.y0,
    (posn.x1-posn.x0)/2, posn.height])
    outline_effect = [patheffects.withStroke(linewidth=3, foreground='black')]
    cbar.set_ticks([])
    cbar.ax.set_xticklabels([])
    #cbar.set_label(f'Speed ({data.variables[u_name].units})',size=18)
    for i in range(10,111,15):
        cbar.ax.text((i/100), .4, i, ha='center', va='center',path_effects=outline_effect,color="w")
    #cbar.ax.text(90, .4, 90, ha='center', va='center',path_effects=outline_effect,color="w")
    
    
    # 500mb Heights
    #---------------------------------------------------------------------------------------------------
    # thredds variables use different isobaric variable names depending on the variable,
    # ie, isobaric1, isobaric7, etc. Thus grab what isobaric number is:
    geo_hgt_iso_num = gfs.get_var_isobaric_num('Geopotential_height_isobaric')
    hgt_500 = hgt[time_strings.index(time),data.variables[geo_hgt_iso_num][:].tolist().index(50000),:,:]
        
    clev500 = np.arange(5200, 6000, 60)
    cs = ax.contour(lons, lats, hgt_500,clev500 ,colors='black', linewidths=2.0,
    linestyles='solid', transform=ccrs.PlateCarree())
    clbls = plt.clabel(cs, **kw_clabels)
    plt.setp(clbls, path_effects=[
        patheffects.withStroke(linewidth=3, foreground="w")])
       
   
    # Save Figure
    #---------------------------------------------------------------------------------------------------    
    WND = im_save_path+"GFS/Wind/"
    if not os.path.isdir(WND):
        os.makedirs(WND)
        
    time_index *= 3
    if time_index < 10:
        times = f"0{time_index}"
    else:
        times = f"{time_index}"
        
    outfile = f"{WND}GFS_0p5_WindSpeed_{file_time}_F{times}.png"
    fig.savefig(outfile,bbox_inches='tight',dpi=120)
    plt.close(fig)