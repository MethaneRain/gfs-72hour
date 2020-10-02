#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 18:47:19 2020

@author: Justin Richling
"""

def datetime_difference(date1,date2):
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