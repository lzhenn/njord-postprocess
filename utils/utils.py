#/usr/bin/env python
"""Commonly used utilities
    Function    
    ---------------
    throw_error(msg):
    write_log(msg, lvl=20):
"""
import sys
import numpy as np
import logging

CWD=sys.path[0]

def throw_error(msg):
    '''
    throw error and exit
    '''
    logging.error(msg)
    exit()

def write_log(msg, lvl=20):
    '''
    write logging log to log file
    level code:
        CRITICAL    50
        ERROR   40
        WARNING 30
        INFO    20
        DEBUG   10
        NOTSET  0
    '''
    logging.log(lvl, msg)

def get_closest_idxy(lat2d, lon2d, lat0, lon0):
    """
        Find the nearest idx, idy in lat2d and lon2d for lat0 and lon0
    """
    dis_lat2d=lat2d-lat0
    dis_lon2d=lon2d-lon0
    dis=abs(dis_lat2d)+abs(dis_lon2d)
    idx=np.argwhere(dis==dis.min())[0].tolist() # x, y position
    return idx[0], idx[1]