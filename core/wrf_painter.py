#!/usr/bin/env python3
'''
Base Class: Painter
    Derived Class: WRFPainter
'''
from . import painter
from utils import utils
from netCDF4 import Dataset

import cmaps
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg') 

import sys, os, subprocess
import wrf

# -------Global Envrionmental Vars--------
print_prefix='lib.painter>>'

CWD=sys.path[0]

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
#SHP_LV3='china_coastline.dbf'
# -------Global Envrionmental Vars--------


class WRFPainter(painter.Painter):
    """WRF Painter"""
    def update(self, cfg):
        '''update WRF specific files'''
        self.wrf_num_dom=cfg['WRF']['num_dom']

    def load_data(self, dom_id):
        '''load WRF files according to domain ID'''
        utils.write_log('get wrfout...')
        fn_stream=subprocess.check_output(
                'ls '+self.arch_root+'/'+self.init_ts+'/wrfout_'+dom_id+'*', 
                shell=True).decode('utf-8')
        
        fn_list=fn_stream.split()
        wrf_list=[Dataset(itm) for itm in fn_list]
        self.wrf_idom=dom_id
        self.wrf_list=wrf_list
        self.wrf_num_file=len(wrf_list)

        self.time_frms=wrf.extract_times(
                wrf_list, timeidx=wrf.ALL_TIMES, 
                method='cat', do_xtime=False)
        
        
        lsmask = wrf.getvar(self.wrf_list[0], 'LANDMASK')
        
        # Get the latitude and longitude points
        self.lats, self.lons = wrf.latlon_coords(lsmask)
        
        utils.write_log('fecthed %3d wrfout files' % len(wrf_list))

    def set_canvas_common(self,var):
        if self.wrf_idom=='d01':
             amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
        proj=wrf.get_cartopy(var)
        
        fig = plt.figure(figsize=[10, 8],frameon=True)
        
        # Set projection and plot the main figure
        ax = fig.add_axes([0.08, 0.01, 0.8, 0.94], projection=proj)

        # Set figure extent
        #ax.set_extent([109, 118, 20, 26],crs=ccrs.PlateCarree())

        # plot shp boundaries
        ax.add_geometries(
                amdn_shp, ccrs.PlateCarree(),
                facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)
        ax.coastlines()

        return ax    
    
    def savefig(self, fig, varname, ifrm):
        '''save fig portal'''
        plt.savefig(
                self.fig_root+'/'+self.init_ts
                +'/'+self.wrf_idom+'.'+varname+'.p%03d.' % ifrm, 
                dpi=120, bbox_inches='tight')
   
    def draw2d_map_t2(self, ifrm):
        '''
        draw 2d spatial 2-m temperature map
        ifrm: ith frame in wrf_list
        '''
        varname='T2'
        title_txt=varname+'@'+str(self.time_frms[ifrm])
        utils.write_log('paint '+title_txt)
        var = wrf.getvar(self.wrf_list, varname, timeidx=ifrm)
        
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(-20,45,66)
        
        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var-273.15),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7)
        
        self.savefig(plt, varname, ifrm)


    def draw2d_map_rh2(self, ifrm):
        '''
        draw 2d spatial 2-m reletive humidity map
        ifrm: ith frame in wrf_list
        '''
        pass
    
    def draw2d_map_wind10(self, ifrm):
        '''
        draw 2d spatial 2-m reletive humidity map
        ifrm: ith frame in wrf_list
        '''
        pass
    
    def draw2d_map_pr3h(self):
        '''
        draw 2d spatial 3-hr total precipitation map
        ifrm: ith frame in wrf_list
        '''
        pass

    def draw2d_map_slp(self):
        '''
        draw 2d spatial sea level pressure map
        ifrm: ith frame in wrf_list
        '''
        pass
