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
        varname='rh2'
        title_txt=varname+'@'+str(self.time_frms[ifrm])
        utils.write_log('paint '+title_txt)
        var = wrf.getvar(self.wrf_list, varname, timeidx=ifrm)
        
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(1,100,50)
        
        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7)
        
        self.savefig(plt, varname, ifrm)
    
    def draw2d_map_wind10(self, ifrm):
        '''
        draw 2d spatial uv vector and speed 
        ifrm: ith frame in wrf_list
        '''
        varname='UV10'
        title_txt=varname+'@'+str(self.time_frms[ifrm])
        utils.write_log('paint '+title_txt)
        u = wrf.getvar(self.wrf_list, 'U10', timeidx=ifrm)
        v = wrf.getvar(self.wrf_list, 'V10', timeidx=ifrm)
        var = u
        var.data = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)
        
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(0,50,50)
        
        shad = plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap)
        
        q_mis=8 # wind vector plotting every q_miss grid
        quv = axe.quiver(wrf.to_np(self.lons[::q_mis,::q_mis]),
                   wrf.to_np(self.lats[::q_mis,::q_mis]),
                   wrf.to_np(u[::q_mis,::q_mis]),wrf.to_np(v[::q_mis,::q_mis]),
                   pivot='mid',units='inches',scale=30,
                   scale_units='inches',color="dimgray",
                   width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

        plt.quiverkey(quv, 0.87, 0.45, 10, r'$10 m/s$', labelpos='N',
                       coordinates='figure')
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(shad, ax=ax, shrink=0.7)
        
        self.savefig(plt, varname, ifrm)
    
    def draw2d_map_pr3h(self):
        '''
        draw 2d spatial 3-hr total precipitation map
        ifrm: ith frame in wrf_list
        '''
        varname='pr3h'
        title_txt=varname+'@'+str(self.time_frms[ifrm])
        utils.write_log('paint '+title_txt)
        var = wrf.getvar(self.wrf_list, 'RAINC', timeidx=ifrm)
        var.data = wrf.to_np(wrf.getvar(self.wrf_list, 'RAINNC', timeidx=ifrm))+
                   wrf.to_np(wrf.getvar(self.wrf_list, 'RAINC', timeidx=ifrm))-
                   wrf.to_np(wrf.getvar(self.wrf_list, 'RAINNC', timeidx=(ifrm-3)))-
                   wrf.to_np(wrf.getvar(self.wrf_list, 'RAINC', timeidx=(ifrm-3)))

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
        plt.colorbar(ax=ax, shrink=0.7)
        
        self.savefig(plt, varname, ifrm)

    def draw2d_map_slp(self):
        '''
        draw 2d spatial sea level pressure map
        ifrm: ith frame in wrf_list
        '''
        varname='slp'#'hPa'
        title_txt=varname+'@'+str(self.time_frms[ifrm])
        utils.write_log('paint '+title_txt)
        var = wrf.getvar(self.wrf_list, varname, timeidx=ifrm)
        
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
    
