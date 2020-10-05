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

_,time_strings,time_final = gfs.get_data_times(data,gfs.u_name)
print(time_strings)


print("\n\nchanging dirs for HiLo script...\n")     
os.chdir("/Users/chowdahead/Documents/gfs-72hour/map-scripts")
print(os.getcwd())
print(os.listdir())
print("importing HiLo script...\n") 
import sys
sys.path.append('/Users/chowdahead/Documents/gfs-72hour/map-scripts')
import jets_250 as jets
print("import done...\ntrying map plot function...\n")
for i in range(0,len(time_strings)): #len(time_strings)+1
    jets.jet250_map(data,i,time_strings,time_final)


import HiLo_map as hilo
print("import done...\ntrying map plot function...\n")
for i in range(0,len(time_strings)): #len(time_strings)+1
    hilo.HiLo_thickness_map(data,i,time_strings,time_final)


