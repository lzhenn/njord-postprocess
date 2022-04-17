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
import matplotlib.dates as mdates
matplotlib.use('Agg') 

import sys, os, subprocess, glob
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
        self.utc_flag=cfg['WRF'].getboolean('wrf_flag')

    def load_metadata(self):
        '''load WRF files according to domain ID'''
        utils.write_log('get meta data from wrfout...')
        pass

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
        self.nrow, self.ncol = (
            self.lats.sizes['south_north'],
            self.lats.sizes['west_east'])
            
        utils.write_log('fecthed %3d wrfout files' % len(wrf_list))
    

    def locate_sta_pos(self, stas):
        '''find station position in idom, irow, jcol'''
        wrf_file=Dataset(self.file_list[0])
        for sta in stas:
            res=(wrf.ll_to_xy(wrf_file,sta.lat, sta.lon))
            irow, icol=res.values[0], res.values[1]
            if irow<0 or irow>(self.nrow-1):
                irow=0
            if icol<0 or icol>(self.nrow-1):
                icol=0
            sta.irow, sta.icol=irow, icol
        return stas

    def set_y1y2_axis(self,sta):
        fig = plt.figure(figsize=(10,8)) 
        ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.35])
        ax1.tick_params(axis='y',labelcolor="red",
            labelsize=SMFONT)
        ax1.tick_params(axis='x',labelsize=SMFONT)
       
        ax2 = ax1.twinx()
        ax2.tick_params(axis='y',labelcolor="blue",
            labelsize=SMFONT)
        
        plt.title("%s (%.2fN, %.2fE)"%(sta.name, sta.lat,
            sta.lon), fontsize=MIDFONT)

        return ax1, ax2

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
    def draw_ts_t2rh2(self, stas):
        '''draw t2m and rh2m thru all stations'''
        utils.write_log('Paint t2rh2 time series')
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        ts1_all=wrf.getvar(
                wrf_list,'T2',timeidx=wrf.ALL_TIMES, 
                method='cat')
        ts2_all=wrf.getvar(
                wrf_list,'rh2',timeidx=wrf.ALL_TIMES, 
                method='cat')
        
        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts1 = ts1_all[:, irow, icol]-273.15
            ts2 = ts2_all[:, irow, icol]
            
            # plot
            ax1, ax2 = self.set_y1y2_axis(sta)
            ax1.plot(self.time_frms, ts1, linewidth=2.0, 
                color="red")
            ax1.set_ylabel(r"T2m ($^\circ$C)",color="red",
                fontsize=SMFONT)
            
            ax2.plot(self.time_frms, ts2, linewidth=2.0,
                color="blue")
            ax2.set_ylabel(r"RH2m (%)",color="blue",
                fontsize=SMFONT)
            ax2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            self.savefig_ts(sta.name+'.ts.t2rh2.png')  

    def draw_ts_pr1h(self, stas):
        '''draw 1hr precipitation thru all stations'''
        utils.write_log('Paint pr1h time series')
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        rainc = wrf.getvar(
                wrf_list,'RAINC',timeidx=wrf.ALL_TIMES, 
                method='cat')
        rainnc= wrf.getvar(
                wrf_list,'RAINNC',timeidx=wrf.ALL_TIMES, 
                method='cat')
        ts1_all = wrf.to_np(rainc[1:,:,:])+wrf.to_np(rainnc[1:,:,:])\
            - wrf.to_np(rainc[0:-1,:,:])-wrf.to_np(rainnc[0:-1,:,:])
        del rainc, rainnc

        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts1 = ts1_all[:, irow, icol]
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            ax1.bar(self.time_frms[1:], ts1, width=0.03,
                edgecolor='black',linewidth=1, align='edge')
            
            ax1.tick_params(axis='both',labelsize=SMFONT)
            ax1.set_ylabel(r"Pr (mm/hr)",fontsize=SMFONT)
            ax1.set_ylim(bottom=0)
            ax1.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("%s (%.2fN, %.2fE)"%(sta.name, sta.lat,
                sta.lon), fontsize=MIDFONT)
            self.savefig_ts(sta.name+'.ts.pr1h.png')  
    
    def draw_ts_wswd10(self, stas):
        '''draw 10-m wind speed and dir thru all stations'''
        utils.write_log('Paint wswd10 time series')
        wdir = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "n"]
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        ts1_all, ts2_all = wrf.getvar(
                wrf_list,'wspd_wdir10',timeidx=wrf.ALL_TIMES, 
                method='cat')
        
        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts1 = ts1_all[:, irow, icol]
            ts2 = ts2_all[:, irow, icol]
            
            # plot
            ax1, ax2 = self.set_y1y2_axis(sta)
            ax1.fill_between(self.time_frms, ts1, color="red")
            ax1.set_ylabel(r"wind 10m (m/s)",color="red",
                fontsize=SMFONT)
            
            ax2.plot(self.time_frms, ts2, "o",
                color="blue")
            ax2.set_yticks(np.arange(0,361,45))
            ax2.set_yticklabels(wdir)
            ax2.set_ylabel(r"wind direction",color="blue",
                fontsize=SMFONT)
            ax2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
        
            self.savefig_ts(sta.name+'.ts.wswd10.png')  
    
    def draw_ts_swslp(self, stas):
        '''draw surface shortwave and slp thru all stations'''
        utils.write_log('Paint swslp time series')
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        ts1_all=wrf.getvar(
                wrf_list,'SWDOWN',timeidx=wrf.ALL_TIMES, 
                method='cat')
        ts2_all=wrf.getvar(
                wrf_list,'slp',timeidx=wrf.ALL_TIMES, 
                method='cat')
        
        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts1 = ts1_all[:, irow, icol]
            ts2 = ts2_all[:, irow, icol]
            
            # plot
            ax1, ax2 = self.set_y1y2_axis(sta)
            ax1.fill_between(self.time_frms, ts1, color="blue")
            ax1.set_ylabel(r"Solar Rad ($W/m^2$)",color="blue",
                fontsize=SMFONT)
            
            ax2.plot(self.time_frms, ts2, "o", 
                color="red")
            ax2.set_ylabel(r"SLP (hPa)",color="blue",
                fontsize=SMFONT)
            ax2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
        
            self.savefig_ts(sta.name+'.ts.swslp.png')  
    
    def draw_hov_cloud(self, stas):
        '''draw hovmoller cloud fraction thru all stations'''
        utils.write_log('Paint hovmoller cloud fraction')
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        term = wrf.getvar(
                wrf_list,'CLDFRA',timeidx=wrf.ALL_TIMES, 
                method='cat')
        pres = wrf.getvar(
                wrf_list,'pressure',timeidx=wrf.ALL_TIMES, 
                method='cat')
        plev = [1000,925,850,700,600, 500,400,350,300,250, 
            200,175,150,125,100]
        dplev = [1000,700, 500,350,250,200,150,100]
        ts_all = wrf.interplevel(term, pres, plev)
        del term, pres
        
        cmap = matplotlib.colors.ListedColormap(
            cmaps.WhiteBlue(np.arange(254,0,-12)))
        levels = np.arange(0,1.01,0.05)
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')
        
        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts = ts_all[:, :, irow, icol]
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            ax = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            
            cont = ax.contourf(self.time_frms,plev, ts.transpose(),
                levels=levels, cmap=cmap, extend='both',norm=norm)
            ax.set_ylabel(r"Pressure (hPa)",fontsize=SMFONT)
            ax.set_yscale("log")
            ax.invert_yaxis()
            ax.set_yticks(dplev,minor=False)
            ax.set_yticklabels(dplev,minor=False)
            plt.minorticks_off()
            
            ax.tick_params(axis='both',labelsize=SMFONT)
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("Cloud Fraction @ %s (%.2fN, %.2fE)"%(
                sta.name, sta.lat, sta.lon), fontsize=MIDFONT)
            plt.colorbar(cont,ax=ax, shrink=0.9,extendfrac='auto')

            self.savefig_ts(sta.name+'.hov_cloud.png')  
    
    def draw_hov_rhwind(self, stas):
        '''draw hovmoller rh and wind thru all stations'''
        utils.write_log('Paint hovmoller rhwind')
        wrf_list=[Dataset(itm) for itm in self.file_list]
        
        term = wrf.getvar(
                wrf_list,'rh',timeidx=wrf.ALL_TIMES, 
                method='cat')
        uterm = wrf.getvar(
                wrf_list,'ua',timeidx=wrf.ALL_TIMES, 
                method='cat')
        vterm = wrf.getvar(
                wrf_list,'va',timeidx=wrf.ALL_TIMES, 
                method='cat')
        pres = wrf.getvar(
                wrf_list,'pressure',timeidx=wrf.ALL_TIMES, 
                method='cat')
        plev = [1000,925,850,700,600, 500,400,350,300,250, 
            200,175,150,125,100]
        dplev = [1000,700, 500,350,250,200,150,100]
        ts_all = wrf.interplevel(term, pres, plev)
        u_all = wrf.interplevel(uterm, pres, plev)
        v_all = wrf.interplevel(vterm, pres, plev)
        del term, pres, uterm, vterm,
        
        cmap = matplotlib.colors.ListedColormap(
            cmaps.MPL_BrBG(np.arange(0,128,6)))
        levels = np.arange(0,101,5)
        norm = matplotlib.colors.BoundaryNorm(boundaries=levels, 
            ncolors=cmap.N,extend='both')
        
        # loop thru stations
        for sta in stas:
            irow, icol=sta.irow, sta.icol
            ts = ts_all[:, :, irow, icol]
            u  =  u_all[:, :, irow, icol]
            v  =  v_all[:, :, irow, icol]
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            ax = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            
            cont = ax.contourf(self.time_frms, plev, ts.transpose(),
                levels=levels, cmap=cmap, extend='both',norm=norm)
            ax.barbs(self.time_frms,plev, 
                u.transpose(), v.transpose(),length=4,
                sizes=dict(emptybarb=0, spacing=0.25, height=0.5),
                barb_increments=dict(half=2, full=4, flag=20))
            ax.set_ylabel(r"Pressure (hPa)",fontsize=SMFONT)
            ax.set_yscale("log")
            ax.invert_yaxis()
            ax.set_yticks(dplev,minor=False)
            ax.set_yticklabels(dplev,minor=False)
            plt.minorticks_off()
            
            ax.tick_params(axis='both',labelsize=SMFONT)
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("RH & wind @ %s (%.2fN, %.2fE)"%(
                sta.name, sta.lat, sta.lon), fontsize=MIDFONT)
            plt.colorbar(cont, ax=ax, shrink=0.9,extendfrac='auto')

            self.savefig_ts(sta.name+'.hov_rhwind.png')  

    def savefig_ts(self, figname):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
        plt.savefig(
                save_dir+'/'+figname, 
                dpi=120, bbox_inches='tight')
        plt.close()

#---------------- Draw 2D Spatial Plot ----------------
    def get_single_var2d(self, varname, ifrm):
        '''get single var 2D in ifrm from wrf file'''
        wrf_file=Dataset(self.file_list[ifrm])
        var = wrf.getvar(wrf_file, varname, timeidx=0)
        return var

    def set_canvas_common(self,var):
        '''set common properties of the canvas'''

        proj=wrf.get_cartopy(var)
        #lats, lons = wrf.latlon_coords(var)
        lats, lons = self.lats, self.lons
        
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
    
    def savefig(self, varname, ifrm):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        fig_prefix=self.wrf_idom+'.'+varname
        
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
       
        plt.savefig(
                save_dir+'/'+fig_prefix+'.p%03d.' % ifrm, 
                dpi=120, bbox_inches='tight')
        plt.close()
    def draw2d_map_sst(self, ifrm, itsk=0):
        '''
        draw 2d spatial SST map
        ifrm: ith frame in wrf_list
        '''
        varname='SST'


        unit='$^\circ$C'
        title_txt=self.wrf_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+self.time_frms[ifrm].strftime('%Y-%m-%d %HZ')
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))
        
        var = self.get_single_var2d(varname, ifrm)
        lsmask = self.get_single_var2d('LANDMASK', ifrm)
        var=var.where(lsmask==0)
        ax=self.set_canvas_common(var)
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(27,31,41)
        
        plt.contourf(
                wrf.to_np(self.lons), wrf.to_np(self.lats), 
                wrf.to_np(var-273.15),
                levels=levels, extend='both', 
                transform=ccrs.PlateCarree(), cmap=cmap)
        
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.7)
        
        self.savefig(varname, ifrm)


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
        
        self.savefig(varname, ifrm)


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
        
        self.savefig(varname, ifrm)
    
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
        
        self.savefig(varname, ifrm)
    
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
        
        self.savefig(varname, ifrm)

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
        
        self.savefig(varname, ifrm)
    
