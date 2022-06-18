#!/usr/bin/env python3
'''
Base Class: Painter
    Derived Class: SWANPainter
'''
from . import painter
from utils import utils
from netCDF4 import Dataset

import cmaps
import numpy as np
import scipy.io as sio
import xarray as xr
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import shapely.geometry as sgeom

import matplotlib
from matplotlib import colors
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
matplotlib.use('Agg') 

import sys, os, subprocess, glob

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


class SWANPainter(painter.Painter):
    """SWAN Painter"""


    def update(self, cfg):
        '''update SWAN specific files'''

        self.swan_num_dom=int(cfg['SWAN']['num_dom'])
        self.grid_temp_dir=cfg['SWAN']['grid_temp_dir']


    def load_metadata(self, dom_id):
        '''load SWAN meta according to domain ID'''

        utils.write_log('get swan grid template file...')
        fn_grid=self.grid_temp_dir+'/roms_'+dom_id+'.nc'
        grd_swan=xr.load_dataset(fn_grid)
        # Get the latitude and longitude points
        self.lats, self.lons = grd_swan['lat_rho'].values, grd_swan['lon_rho'].values
        self.nrow, self.ncol = self.lats.shape

        # land sea mask
        self.lsmask =grd_swan['mask_rho']
        
        utils.write_log('get swan hist file...')
        fn_stream=subprocess.check_output(
            'ls '+self.arch_root+'/'+self.init_ts+'/hsig_*.mat', 
            #'ls '+self.arch_root+'/'+self.init_ts+'/hsig_*_'+dom_id+'.mat', 
            shell=True).decode('utf-8')
        
        fn_list=fn_stream.split()
        self.file_list=fn_list
        self.swan_idom=dom_id
        self.swan_num_file=len(fn_list)
    

    def load_data(self, fid):
        '''load SWAN data according to fid and domain ID'''
        
        mat_swan=sio.loadmat(self.file_list[fid])
        self.swan_mat_key=[]
        for key in mat_swan:
            if 'Hsig' in key:
                self.swan_mat_key.append(key)
        self.mat_swan=mat_swan
    
    def load_ts_data(self, varname, sta):
        '''load SWAN data according to fid and domain ID'''
        
        idx, idy=sta.irow, sta.icol
        self.ts1_all=[]
        self.time_frms=[]
        for fid in range(self.swan_num_file):
            mat_swan=sio.loadmat(self.file_list[fid])
            for key in mat_swan:
                if varname in key:
                    trfm_lst=key.split('_')
                    ymdhms=trfm_lst[1]+'_'+trfm_lst[2]
                    tfrm_dt=datetime.strptime( ymdhms,'%Y%m%d_%H%M%S')
                    if tfrm_dt not in self.time_frms:
                        self.time_frms.append(tfrm_dt)
                        self.ts1_all.append(mat_swan[key][idx, idy])
    
    def locate_sta_pos(self, stas):
        '''find station position in irow, icol, at sea or not'''
        lat2d=self.lats
        lon2d=self.lons
        for sta in stas:
            irow, icol=utils.get_closest_idxy(
                lat2d, lon2d, sta.lat, sta.lon)
            if irow<0 or irow>(self.nrow-1):
                irow=0
            if icol<0 or icol>(self.ncol-1):
                icol=0
            sta.irow, sta.icol=irow, icol
            sta.is_sea=self.lsmask[irow, icol]
        return stas

    def form_anim(self, opt=''):
        '''form animation from pngs'''
        save_dir=self.fig_root+'/'+self.init_ts
        prefix2d_lst=[]
        lst_png=glob.glob(save_dir+'/*.png')
        for png in lst_png:
            fn=png.split('/')[-1]
            fn=fn.split('.')[0]+'.'+fn.split('.')[1]
            if not(fn in prefix2d_lst):
                prefix2d_lst.append(fn)
        if opt in prefix2d_lst:
            prefix2d_lst=[opt]
        
        for itm in prefix2d_lst:
            utils.write_log('form animation for %s' % itm)
            cmd='convert -delay 30 -loop 0 '+save_dir+'/'+itm+'*.png '+save_dir+'/'+itm+'.gif'
            subprocess.call(cmd, shell=True) 


 #---------------- Draw Station Time Series ----------------
    def draw_ts_hsig(self, stas):
        '''draw sigwave height and peak wave period thru all stations'''
        # loop thru stations
        varname='Hsig'
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint '+varname+' time series for %s' % sta.name)
            self.load_ts_data(varname, sta)
            # plot
            fig = plt.figure(figsize=(10,8)) 
            ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            ax1.plot(self.time_frms, self.ts1_all, linewidth=2.0) 
            
            ax1.tick_params(axis='both',labelsize=SMFONT)
            ax1.set_ylabel(r"Significant Wave Height (m)",fontsize=SMFONT)
            ax1.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("%s (%.2fN, %.2fE)"%(sta.name, sta.lat,
                sta.lon), fontsize=MIDFONT)
            self.savefig_ts(sta.name+'.ts.'+varname+'.png')  


    def draw_ts_zeta(self, stas):
        '''draw free sea surf height thru all stations'''
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint zeta time series for %s' % sta.name)
            ts1_all=self.ds_swan['zeta'][:,sta.irow,sta.icol]
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            ax1.plot(self.time_frms, ts1_all, linewidth=2.0) 
            
            ax1.tick_params(axis='both',labelsize=SMFONT)
            ax1.set_ylabel(r"sea surf height (m)",fontsize=SMFONT)
            ax1.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("%s (%.2fN, %.2fE)"%(sta.name, sta.lat,
                sta.lon), fontsize=MIDFONT)
            self.savefig_ts(sta.name+'.ts.zeta.png')  


#---------------- Draw 2D Spatial Plot ----------------
    def draw_2d_map_hsig(self, trfm, itsk=0):
        '''
        draw 2d spatial significant wave height map
        '''
        varname='Hsig'
        unit='m'
        var2d=self.mat_swan[trfm]

        trfm_lst=trfm.split('_')
        ymdhms=trfm_lst[1]+'_'+trfm_lst[2]
        tfrm_dt=datetime.strptime( ymdhms,'%Y%m%d_%H%M%S')
        tfrm_str=tfrm_dt.strftime('%Y-%m-%d %H:%M:%S')
        title_txt=self.swan_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        cnlevels=np.linspace(0,7.,71)

        ncmap = ListedColormap(cmap(range(0,200,1))) #cmaps.precip2_17lev
        #ncmap = ListedColormap(["red","blue","green"])
        #cnlevels = np.concatenate((np.arange(-63,0,1),np.arange(0,615,15)))
        norm1  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
        '''
        plt.contourf(
            self.lons, self.lats, 
            var2d,
            levels=levels, extend='both', norm=norm1, 
            transform=ccrs.PlateCarree(), cmap=cmap)
        '''
        plt.pcolormesh(
            self.lons, self.lats, 
            var2d,
            norm=norm1, 
            transform=ccrs.PlateCarree(), cmap=cmap)

        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm_dt.strftime('%Y%m%d%H'))

    def set_canvas_common(self):
        '''set common properties of the canvas'''
        lonmin, lonmax=self.lons.min(), self.lons.max()
        latmin, latmax=self.lats.min(), self.lats.max()

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
        
 
        amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV4).geometries()
        amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
        #ax.add_geometries(
        #    amdn_shp_outer, ccrs.PlateCarree(),
        #    facecolor='none', edgecolor='black',linewidth=1, zorder = 1)
        #ax.coastlines()
        
        # plot shp boundaries
        ax.add_geometries(
            amdn_shp, ccrs.PlateCarree(),
            facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)
    
        return ax    
    
    def savefig(self, varname, time_str):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        fig_prefix='swan.'+self.swan_idom+'.'+varname
        
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
        plt.savefig(
                save_dir+'/'+fig_prefix+'.p%s.' % time_str, 
                dpi=120, bbox_inches='tight')
        plt.close()

