#!/usr/bin/env python3
'''
ROMS UV 

History:
Apr 16, 2022 --- Zhenning LI

'''
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
matplotlib.use('Agg') 


# Constants
BIGFONT=22
MIDFONT=18
SMFONT=14

def get_closest_idxy(lat2d, lon2d, lat0, lon0):
    '''
        Find the nearest idx, idy in lat2d and lon2d for lat0 and lon0
    ''' 
    dis_lat2d=lat2d-lat0
    dis_lon2d=lon2d-lon0
    dis=abs(dis_lat2d)+abs(dis_lon2d)
    idx=np.argwhere(dis==dis.min())[0].tolist() # x, y position
    return idx[0], idx[1]

if __name__=='__main__':
    hist_file='/home/metctm1/array/data/1911-COAWST/ERA5_TY2001_org/gba_ocean_his.nc'
    
    ds_all=xr.open_dataset(hist_file)
    var=ds_all['zeta']
    var2=ds_all['Hwave']
    lat0=22.2375
    lon0=114.30
    idx, idy=get_closest_idxy(ds_all.lat_rho.values, ds_all.lon_rho.values, lat0, lon0)
    istrt=48
    iend=97
    fig = plt.figure(figsize=(10,8)) 
    ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.35])
    ax1.plot(var.ocean_time[istrt:iend], 
        0.2*var2[istrt:iend,idx,idy]+var[istrt:iend,idx,idy], 
        linewidth=2.0, color='k', label='Total Surge') 
    ax1.plot(var.ocean_time[istrt:iend], 0.2*var2[istrt:iend,idx,idy], 
        linewidth=2.0, color='blue', label='Wave Setup') 
    ax1.legend(loc='upper right', fontsize=SMFONT) 
    ax1.tick_params(axis='both',labelsize=SMFONT)
    ax1.set_ylabel(r"Sea level surge (m)",fontsize=SMFONT)
    ax1.xaxis.set_major_formatter(
        mdates.DateFormatter("%b %d\n%H:00"))
    
    plt.title("%s (%.2fN, %.2fE)"%('Querry Bay', lat0,
        lon0), fontsize=MIDFONT)
    plt.savefig(
                'zeta.png', 
                dpi=120, bbox_inches='tight')
    plt.close()
    ds_all.close()
