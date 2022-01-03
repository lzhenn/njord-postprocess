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
import pandas as pd
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
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
        self.file_list=fn_list
        wrf_list=[Dataset(itm) for itm in fn_list]
        self.wrf_idom=dom_id
        self.wrf_num_file=len(wrf_list)

        self.time_frms=wrf.extract_times(
                wrf_list, timeidx=wrf.ALL_TIMES, 
                method='cat', do_xtime=False)
        self.time_frms = pd.to_datetime(
            np.datetime_as_string(self.time_frms,unit='h'),
            format='%Y-%m-%dT%H')
        
        lsmask = wrf.getvar(wrf_list[0], 'LANDMASK')
        
        # Get the latitude and longitude points
        self.lats, self.lons = wrf.latlon_coords(lsmask)
        
        utils.write_log('fecthed %3d wrfout files' % len(wrf_list))

    def set_canvas_common(self,var):
        if self.wrf_idom=='d01':
             amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
             lat_sp = 5.0 
             lon_sp = 5.0

        proj=wrf.get_cartopy(var)
        lats, lons = wrf.latlon_coords(var)
        
        fig = plt.figure(figsize=[10, 8],frameon=True)
        
        # Set projection and plot the main figure
        ax = fig.add_axes([0.08, 0.01, 0.8, 0.94], projection=proj)

        # Set figure extent
        #ax.set_extent([109, 118, 20, 26],crs=ccrs.PlateCarree())
        #ax.set_xticks(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),lat_sp), crs=ccrs.PlateCarree())
        #ax.set_yticks(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),lon_sp), crs=ccrs.PlateCarree())
        #ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
        #ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        #ax.xaxis.set_major_formatter(LongitudeFormatter())
        #ax.yaxis.set_major_formatter(LatitudeFormatter())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                color='grey', linestyle='--', linewidth=1)
        gl.xlabels_top = False
        gl.ylabels_right = False
        #gl.xlines = False
        #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
        #gl.xformatter = LONGITUDE_FORMATTER
        #gl.yformatter = LATITUDE_FORMATTER
        #gl.xlabel_style = {'size': 15, 'color': 'gray'}
        #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

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

    def get_single_var2d(self, varname, ifrm):
        '''get single var 2D in ifrm from wrf file'''
        wrf_file=Dataset(self.file_list[ifrm])
        var = wrf.getvar(wrf_file, varname, timeidx=0)
        return var

    def draw2d_map_t2(self, ifrm, itsk=0):
        '''
        draw 2d spatial 2-m temperature map
        ifrm: ith frame in wrf_list
        '''
        varname='T2'
        title_txt=varname+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        
        var = self.get_single_var2d(varname, ifrm)
        
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


    def draw2d_map_rh2(self, ifrm, itsk=0):
        '''
        draw 2d spatial 2-m reletive humidity map
        ifrm: ith frame in wrf_list
        '''
        varname='rh2'
        title_txt=varname+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        
        var = self.get_single_var2d(varname, ifrm)
        
        ax=self.set_canvas_common(var)
        
        cmap=matplotlib.colors.ListedColormap(
            cmaps.CBR_drywet(
            np.concatenate((np.arange(0,5,1),np.arange(6,11,1)))))
        levels=np.arange(20,110,10)
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')
        
        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7,extendfrac='auto')
        
        self.savefig(plt, varname, ifrm)
    
    def draw2d_map_wind10(self, ifrm, itsk=0):
        '''
        draw 2d spatial uv vector and speed 
        ifrm: ith frame in wrf_list
        '''
        varname='UV10'
        title_txt=varname+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        u = self.get_single_var2d('U10', ifrm)
        v = self.get_single_var2d('V10', ifrm)
        
        var = self.get_single_var2d('slp', ifrm)
        var.data = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)
        
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.wind_17lev
        levels=np.arange(0,34,2)
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')
        
        shad = plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
        
        q_mis=12 # wind vector plotting every q_miss grid
        quv = ax.quiver(wrf.to_np(self.lons[::q_mis,::q_mis]),
                   wrf.to_np(self.lats[::q_mis,::q_mis]),
                   wrf.to_np(u[::q_mis,::q_mis]),wrf.to_np(v[::q_mis,::q_mis]),
                   pivot='mid',units='inches',scale=40,
                   scale_units='inches',color="dimgray",
                   width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

        plt.quiverkey(quv, 0.77, 0.1, 10, r'$10 m/s$', labelpos='N',
                       coordinates='figure')
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(shad, ax=ax, shrink=0.7,extendfrac='auto')
        
        self.savefig(plt, varname, ifrm)
    
    def draw2d_map_pr3h(self, ifrm, itsk=0):
        '''
        draw 2d spatial 3-hr total precipitation map
        ifrm: ith frame in wrf_list
        '''
        if ifrm<3:
            return
        varname='pr3h'
        title_txt=varname+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        var0c = self.get_single_var2d('RAINC', ifrm)
        var0nc = self.get_single_var2d('RAINNC', ifrm)
        var3c = self.get_single_var2d('RAINC', ifrm-3)
        var3nc = self.get_single_var2d('RAINNC', ifrm-3)
        var = var0c
        var.data = (wrf.to_np(var0c) + wrf.to_np(var0nc) - wrf.to_np(var3c)
                - wrf.to_np(var3nc))

        ax=self.set_canvas_common(var)
        
        cmap = cmaps.precip2_17lev
        levels = [0.1,0.5,1,3,5,10,15,20,30,40,60,80,100,120,150,200,250] 
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')

        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7,extendfrac='auto')
        
        self.savefig(plt, varname, ifrm)

    def draw2d_map_slp(self, ifrm, itsk=0):
        '''
        draw 2d spatial sea level pressure map
        ifrm: ith frame in wrf_list
        '''
        varname='slp'#'hPa'
        title_txt=varname+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        var = self.get_single_var2d(varname, ifrm)
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(950,1050,101)
        
        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7)
        
        self.savefig(plt, varname, ifrm)
    
