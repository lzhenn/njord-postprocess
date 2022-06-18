#!/usr/bin/env python3
'''
ROMS zeta interpolator 

History:
Apr 24, 2022 --- Zhenning LI

'''
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata
import sys

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import cmaps
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
matplotlib.use('Agg') 


# Constants
BIGFONT=22
MIDFONT=18
SMFONT=14
# Province level
SHP_LV1='cnhimap.dbf'
# City level
SHP_LV2='gadm36_CHN_2.dbf'
# County level
SHP_LV3='county_2004.dbf'
# For D03 Coastline
SHP_LV4='china_coastline.dbf'


def set_canvas_common(grd):
    '''set common properties of the canvas'''
    lonmin, lonmax=grd['lon_rho'].min(), grd['lon_rho'].max()
    latmin, latmax=grd['lat_rho'].min(), grd['lat_rho'].max()

    proj = ccrs.Mercator(
        central_longitude=115., min_latitude=-80.0, max_latitude=84.0, 
        globe=None, latitude_true_scale=23.0, false_easting=0.0, 
        false_northing=0.0, scale_factor=None) 
    
    fig = plt.figure(figsize=[10, 8],frameon=True)
    
    # Set projection and plot the main figure
    ax = fig.add_axes([0.08, 0.01, 0.8, 0.94], projection=proj)

    # Set figure extent
    ax.set_extent(
        [lonmin, lonmax, latmin, latmax],
        crs=ccrs.PlateCarree())

    # add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
        color='grey', linestyle='--', linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False
    

    CWD=sys.path[0]+'/../'
    amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV4).geometries()
    amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
    ax.add_geometries(
        amdn_shp_outer, ccrs.PlateCarree(),
        facecolor='none', edgecolor='black',linewidth=1, zorder = 1)
    #ax.coastlines()
    
    # plot shp boundaries
    ax.add_geometries(
        amdn_shp, ccrs.PlateCarree(),
        facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)

    return ax    


def savefig(varname, time_str):
    '''save fig portal'''
    save_dir='./'
    fig_prefix='swan.'+varname
    
    plt.savefig(
            save_dir+'/'+fig_prefix+'.p%s.' % time_str, 
            dpi=120, bbox_inches='tight')
    plt.close()


if __name__=='__main__':
    tgt_dom_file='/home/lzhenn/array74/workspace/calypso_pipeline/domaindb/swant1t2/roms_d02.nc'
    src_file='/home/metctm1/array/data/1911-COAWST/ERA5_TY2001_org/gba_ocean_his.nc'
    
    varname='surge_height'
    unit='m'
    grd_swan=xr.load_dataset(tgt_dom_file)
    # Get the latitude and longitude points
    tgt_lats, tgt_lons = grd_swan['lat_rho'].values, grd_swan['lon_rho'].values
    tgt_nrow, tgt_ncol = tgt_lats.shape

    grd_roms=xr.open_dataset(src_file)
    # Get the latitude and longitude points
    src_lats, src_lons = grd_roms['lat_rho'].values, grd_roms['lon_rho'].values
    src_nrow, src_ncol = src_lats.shape
    
    for test_frm in range(48,145,2):
        var=grd_roms['zeta'][test_frm,:,:]
        tfrm=grd_roms['ocean_time'][test_frm]
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        points=src_lats.ravel(), src_lons.ravel()
        values=var.values.ravel()
        grid_z0 = griddata(
            points, values, (tgt_lats.ravel(), tgt_lons.ravel()), 
            method='linear')
        df = pd.DataFrame(grid_z0 )
        grid_z0=df.interpolate(method='linear').values
        tgt_var = grid_z0.reshape(tgt_nrow, tgt_ncol)
        tgt_var[grd_swan['mask_rho']==0]=np.nan
        ax=set_canvas_common(grd_swan)
        
        cmap=cmaps.cmocean_deep
        #cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(0,5.,101)
        plt.contourf(
            tgt_lons, tgt_lats, 
            tgt_var,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        
        title_txt='d02: '+varname+' ('+unit+') @'+tfrm_str
        print(title_txt)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    grd_roms.close() 



