#!/usr/bin/env python3
'''
ROMS Vertical Profile Utilities

History:
Apr 16, 2022 --- Zhenning LI

'''
import numpy as np
import xarray as xr

def get_closest_idxy(lat2d, lon2d, lat0, lon0):
    '''
        Find the nearest idx, idy in lat2d and lon2d for lat0 and lon0
    ''' 
    dis_lat2d=lat2d-lat0
    dis_lon2d=lon2d-lon0
    dis=abs(dis_lat2d)+abs(dis_lon2d)
    idx=np.argwhere(dis==dis.min())[0].tolist() # x, y position
    return idx[0], idx[1]

def get_ds_with_mpt_dpth(his_file, lat, lon):
    '''
        Get midpoint depths for ROMS output 
        return in dataset
    '''
    ds = xr.open_dataset(his_file)

    # vertical transform
    if ds.Vtransform == 1:
        Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho = Zo_rho + ds.zeta * (1 + Zo_rho / ds.h)
    elif ds.Vtransform == 2:
        Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
    # z_rho [layer, x, y, time]
    ds.coords["z_rho"] = z_rho.transpose()  # needing transpose seems to be an xarray bug
    return ds
 

def integ_uvbar(ds, idx, idy, max_depth=-1):
    '''
        Integrate ubar and vbar to get ubar_int and vbar_int
    '''
    
    # z_rho [time, layer]
    z_rho=ds.z_rho[:,idx,idy,:].values.transpose()
    
    # get depth on each model layer
    dh=z_rho[:,1:]-z_rho[:,:-1]
    
    # get ubar and vbar [time, layer, y, x]
    ubar,vbar=ds.u[:,:, idy, idx], ds.v[:,:, idy, idx]
    
    # mask uv deeper than max depth
    ubar=xr.where(z_rho<=max_depth, 0, ubar)
    vbar=xr.where(z_rho<=max_depth, 0, vbar)
    # integrate ubar and vbar
    ubar_int=np.sum(ubar[:,1:]*dh, axis=1)
    vbar_int=np.sum(vbar[:,1:]*dh, axis=1)
    print(ubar_int)
    print(vbar_int)
    return ubar_int, vbar_int

if __name__=='__main__':
    hist_file='/home/lzhenn/cooperate/data/case_study/coupled/2021091512/njord_his_d01.20210915.nc'
    lat=20.961667
    lon=111.608889    
    
    # add depth (m) to dataset
    ds=get_ds_with_mpt_dpth(hist_file,lat,lon)

    # get closest idx, idy
    idx, idy = get_closest_idxy(
        ds.lat_rho.values, ds.lon_rho.values, lat, lon)

    # integrate ubar and vbar
    max_depth=-5
    ubar_int, vbar_int=integ_uvbar(ds, idx, idy, max_depth)