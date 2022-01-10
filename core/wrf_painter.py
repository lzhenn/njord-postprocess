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
# For D03 Coastline
SHP_LV4='china_coastline.dbf'
# -------Global Envrionmental Vars--------


class WRFPainter(painter.Painter):
    """WRF Painter"""
    def update(self, cfg):
        '''update WRF specific files'''
        self.wrf_num_dom=int(cfg['WRF']['num_dom'])

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
        '''set common properties of the canvas'''

        proj=wrf.get_cartopy(var)
        lats, lons = wrf.latlon_coords(var)
        
        fig = plt.figure(figsize=[10, 8],frameon=True)
        
        # Set projection and plot the main figure
        ax = fig.add_axes([0.08, 0.01, 0.8, 0.94], projection=proj)

        # add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                color='grey', linestyle='--', linewidth=1)
        gl.top_labels = False
        gl.right_labels = False
        
 
        if self.wrf_idom=='d01':
            amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
            ax.coastlines()
        elif self.wrf_idom=='d02':
            amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV3).geometries()
            amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
            ax.add_geometries(
                    amdn_shp_outer, ccrs.PlateCarree(),
                    facecolor='none', edgecolor='black',linewidth=1, zorder = 1)
            ax.coastlines()
        elif self.wrf_idom=='d03':
            amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV3).geometries()
            amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV4).geometries()
            ax.add_geometries(
                    amdn_shp_outer, ccrs.PlateCarree(),
                    facecolor='none', edgecolor='black',linewidth=1, zorder = 1)

        # plot shp boundaries
        ax.add_geometries(
                amdn_shp, ccrs.PlateCarree(),
                facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)
        
        return ax    
    
    def savefig(self, fig, varname, ifrm):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
        plt.savefig(
                save_dir+'/'+self.wrf_idom+'.'+varname+'.p%03d.' % ifrm, 
                dpi=120, bbox_inches='tight')

    def get_single_var2d(self, varname, ifrm):
        '''get single var 2D in ifrm from wrf file'''
        wrf_file=Dataset(self.file_list[ifrm])
        var = wrf.getvar(wrf_file, varname, timeidx=0)
        return var


    def draw_ts_t2rh2(self, loc):
        pass

    def draw2d_map_t2(self, ifrm, itsk=0):
        '''
        draw 2d spatial 2-m temperature map
        ifrm: ith frame in wrf_list
        '''
        varname='T2'


        unit='$^\circ$C'
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        
        var = self.get_single_var2d(varname, ifrm)
        
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(-10,40,101)
        
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
        unit='%'
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        
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
        
        unit='m/s'
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        u = self.get_single_var2d('U10', ifrm)
        v = self.get_single_var2d('V10', ifrm)
        
        wspd10 = self.get_single_var2d('wspd_wdir10', ifrm)[0,:,:]
        ax=self.set_canvas_common(wspd10)
        
        cmap=cmaps.wind_17lev
        levels=np.arange(0,34,2)
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')
        
        shad = plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(wspd10),
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
    
    def draw2d_map_pr(self, ifrm, hrs, itsk=0):
        '''
        draw 2d spatial n-hr total precipitation map
        ifrm: ith frame in wrf_list
        '''
        if ifrm<hrs:
            return
        varname='pr%dh' % hrs
        unit='mm/%dhr' % hrs
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        rainc0 = self.get_single_var2d('RAINC', ifrm)
        rainnc0 = self.get_single_var2d('RAINNC', ifrm)
        rainc_hrs = self.get_single_var2d('RAINC', ifrm-hrs)
        rainnc_hrs = self.get_single_var2d('RAINNC', ifrm-hrs)
        var = rainc0
        var.data = (wrf.to_np(rainc0) + wrf.to_np(rainnc0) - wrf.to_np(rainc_hrs)
                - wrf.to_np(rainnc_hrs))

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
        unit='hPa'
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
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
    
