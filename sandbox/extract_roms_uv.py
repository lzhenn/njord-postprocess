#!/usr/bin/env python3
'''
ROMS UV 

History:
Apr 16, 2022 --- Zhenning LI

'''
import numpy as np
import xarray as xr



if __name__=='__main__':
    hist_file='/home/lzhenn/cooperate/data/case_study/coupled/2021091512/njord_his_d01.20210915.nc'
    
    ds_all=xr.open_dataset(hist_file)
    ds_all['u'].to_netcdf('/home/lzhenn/cooperate/data/case_study/coupled/2021091512/njord_his_d01.u.20210915.nc')
    ds_all.close()
