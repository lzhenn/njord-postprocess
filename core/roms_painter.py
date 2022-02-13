#!/usr/bin/env python3
'''
Base Class: Painter
    Derived Class: ROMSPainter
'''
from . import painter
from utils import utils
from netCDF4 import Dataset

import cmaps
import numpy as np
import pandas as pd
import xarray as xr
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


class ROMSPainter(painter.Painter):
    """ROMS Painter"""
    def update(self, cfg):
        '''update ROMS specific files'''
        self.utc_flag=cfg['ROMS'].getboolean('roms_flag')

    def load_data(self, dom_id):
        '''load ROMS files according to domain ID'''
        utils.write_log('get roms hist file...')
        fn_stream=subprocess.check_output(
            'ls '+self.arch_root+'/'+self.init_ts+'/njord_his_'+dom_id+'*', 
            shell=True).decode('utf-8')
        
        fn_list=fn_stream.split()
        self.file_list=fn_list
        self.roms_idom=dom_id
        self.roms_num_file=len(fn_list)
        ds_roms=xr.open_mfdataset(fn_list)
        # land sea mask
        self.ds_roms=ds_roms
        self.lsmask =ds_roms['mask_rho'][0,:,:]
        # vertical coordinate
        s_rho=ds_roms['s_rho'][:]
        self.isurf=s_rho.size-1
        # Get the latitude and longitude points
        self.lats, self.lons = ds_roms['lat_rho'].values, ds_roms['lon_rho'].values
        self.nrow, self.ncol = self.lats.shape
        self.time_frms=ds_roms['ocean_time'][:]
        self.nfrms_file=(self.time_frms.size-self.roms_num_file)//self.roms_num_file+1
        utils.write_log('fecthed %3d romsout files' % len(fn_list))
        ds_roms.close()
    def locate_sta_pos(self, stas):
        '''find station position in irow, icol, at sea or not'''
        lat2d=self.lats
        lon2d=self.lons
        for sta in stas:
            irow, icol=utils.get_closest_idxy(
                lat2d, lon2d, sta.lat, sta.lon)
            if irow<0 or irow>(self.nrow-1):
                irow=0
            if icol<0 or icol>(self.nrow-1):
                icol=0
            sta.irow, sta.icol=irow, icol
            sta.is_sea=self.lsmask[irow, icol]
        return stas

    def form_anim(self):
        '''form animation from pngs'''
        save_dir=self.fig_root+'/'+self.init_ts
        prefix2d_lst=[]
        lst_png=glob.glob(save_dir+'/*.png')
        for png in lst_png:
            fn=png.split('/')[-1]
            fn=fn.split('.')[0]+'.'+fn.split('.')[1]
            if not(fn in prefix2d_lst):
                prefix2d_lst.append(fn)
        for itm in prefix2d_lst:
            utils.write_log('form animation for %s' % itm)
            cmd='convert -delay 30 -loop 0 '+save_dir+'/'+itm+'*.png '+save_dir+'/'+itm+'.gif'
            subprocess.call(cmd, shell=True) 

 #---------------- Draw Station Time Series ----------------
    def draw_ts_sstsss(self, stas):
        '''draw sea surf temp and sea surf sanility thru all stations'''
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint SST&SSS time series for %s' % sta.name)
            ts1_all=self.ds_roms['temp'][:,self.isurf,sta.irow,sta.icol]
            ts2_all=self.ds_roms['salt'][:,self.isurf,sta.irow,sta.icol]
            
            # plot
            ax1, ax2 = self.set_y1y2_axis(sta)
            ax1.plot(self.time_frms, ts1_all, linewidth=2.0, 
                color="red")
            ax1.set_ylabel(r"SST ($^\circ$C)",color="red",
                fontsize=SMFONT)
            ax2.plot(self.time_frms, ts2_all, linewidth=2.0,
                color="blue")
            ax2.set_ylabel(r"SSS (g/kg)",color="blue",
                fontsize=SMFONT)
            ax2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            self.savefig_ts(sta.name+'.ts.sstsss.png')  

    def draw_ts_hwave_lwavep(self, stas):
        '''draw sigwave height and peak wave period thru all stations'''
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint Hwave&Pwave time series for %s' % sta.name)
            ts1_all=self.ds_roms['Hwave'][:,sta.irow,sta.icol]
            ts2_all=self.ds_roms['Pwave_top'][:,sta.irow,sta.icol]
            
            # plot
            ax1, ax2 = self.set_y1y2_axis(sta)
            ax1.plot(self.time_frms, ts1_all, linewidth=2.0, 
                color="red")
            ax1.set_ylabel(r"Hwave (m)",color="red",
                fontsize=SMFONT)
            ax2.plot(self.time_frms, ts2_all, linewidth=2.0,
                color="blue")
            ax2.set_ylabel(r"Pwave (s)",color="blue",
                fontsize=SMFONT)
            ax2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            self.savefig_ts(sta.name+'.ts.hwavelwavep.png')  

    def draw_ts_zeta(self, stas):
        '''draw free sea surf height thru all stations'''
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint zeta time series for %s' % sta.name)
            ts1_all=self.ds_roms['zeta'][:,sta.irow,sta.icol]
            
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

#---------------- Draw 2D Hovmoller Plot ----------------
    def get_var_timedepth(self, varname, sta, mdepth, sdepth):
        '''get var 2D plane time versus depth 
            varname: variable name
            sta: stations
            mdepth: max depth
            sdepth: depth spacing
        '''
        zeta = self.ds_roms['zeta'][:,sta.irow,sta.icol] 
        ds = xr.open_dataset(self.file_list[0])
        h = ds['h'][sta.irow,sta.icol] 
        if ds.Vtransform == 1:
            Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * h 
            z_rho = Zo_rho + zeta * (1 + Zo_rho / h)
        elif ds.Vtransform == 2:
            Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * h) / (ds.hc + h)
            z_rho = zeta + (zeta + h) * Zo_rho
        ds.close()

        depth = np.arange(mdepth,0,sdepth)
        data = np.empty((len(self.time_frms),len(depth)), dtype = float)
        for nt in range(len(self.time_frms)):
            term = self.ds_roms[varname][nt,:,sta.irow,sta.icol] 
            term.coords['s_rho'] = z_rho[nt,:].values
            data[nt,:] = term.interp(
                s_rho=depth, method="cubic").values

        var = xr.DataArray(data, coords=[
            ("ocean_time", self.time_frms), ("z_rho", depth)])
        return var

    def draw_hov_swt(self, stas, max_depth=-200, sdepth=10):
        '''draw hovmoller sea water temp thru all stations
            stas: list of stations
            max_depth: max depth to plot
            sdepth: depth spacing used to interp
        '''
        cnlevels = np.arange(20,29.5,0.25)
        ncmap = matplotlib.colors.ListedColormap(
            cmaps.MPL_jet(range(0,127,3))) #cmaps.precip2_17lev
        norm = matplotlib.colors.BoundaryNorm(
            boundaries=cnlevels,ncolors=ncmap.N,extend='both')
        
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint hovmoller sea water temp for %s' % sta.name)
            var = self.get_var_timedepth('temp', sta, max_depth, sdepth)
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            axe = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            cont = axe.contourf(var.ocean_time, var.z_rho, var.transpose(), 
                cnlevels,cmap=ncmap,extend='both',norm=norm)
            plt.colorbar(cont, ax=axe)

            axe.set_ylim([max_depth,0])
            axe.set_ylabel("depth (m)",fontsize=SMFONT)  # Add an x-label to the axes.
            axe.tick_params(axis='both',labelsize=SMFONT)
            axe.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("%s (%.2fN, %.2fE) oceam temp ($^\circ$C)"%(
                sta.name, sta.lat, sta.lon), fontsize=MIDFONT)
            self.savefig_ts(sta.name+'.hov.swt.png')  

    def draw_hov_sws(self, stas, max_depth=-200,sdepth=10):
        '''draw hovmoller sea water salt thru all stations
            stas: list of stations
            max_depth: max depth to plot
            sdepth: depth spacing used to interp
        '''
        cnlevels = np.arange(20,39,0.5) #500Z, gpm
        ncmap = matplotlib.colors.ListedColormap(
            cmaps.MPL_jet(range(0,127,3))) #cmaps.precip2_17lev
        norm = matplotlib.colors.BoundaryNorm(
            boundaries=cnlevels,ncolors=ncmap.N,extend='both')
        
        # loop thru stations
        for sta in stas:
            if sta.is_sea==0:
                utils.write_log('Station name: %s is not at sea' % sta.name)
                continue
            utils.write_log('Paint hovmoller sea water salt for %s' % sta.name)
            var = self.get_var_timedepth('salt',sta, max_depth, sdepth)
            
            # plot
            fig = plt.figure(figsize=(10,8)) 
            axe = fig.add_axes([0.1, 0.1, 0.9, 0.35])
            cont = axe.contourf(var.ocean_time, var.z_rho, var.transpose(), 
                cnlevels,cmap=ncmap,extend='both',norm=norm)
            plt.colorbar(cont, ax=axe)

            axe.set_ylim([max_depth,0])
            axe.set_ylabel("depth (m)",fontsize=SMFONT)  # Add an x-label to the axes.
            axe.tick_params(axis='both',labelsize=SMFONT)
            axe.xaxis.set_major_formatter(
                mdates.DateFormatter("%b %d\n%H:00"))
            
            plt.title("%s (%.2fN, %.2fE) ocean sanility (g/kg)"%(
                sta.name, sta.lat, sta.lon), fontsize=MIDFONT)
            self.savefig_ts(sta.name+'.hov.sws.png')  

#---------------- Draw 2D Spatial Plot ----------------
    def draw2d_map_sst(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial sea surface temp map
        '''
        varname='SST'
        unit='$^\circ$C'

        # read sst
        var2d=self.get_single_var2din3d(fn, 'temp', ifrm, self.isurf)
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(10,35,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    def draw2d_map_sss(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial sea surface salinity map
        '''
        varname='SSS'
        unit='g/kg'

        # read sst
        var2d=self.get_single_var2din3d(fn, 'salt', ifrm, self.isurf)
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(20,45,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    def draw2d_map_hwave(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial sigwave height map
        '''
        varname='hwave'
        unit='m'

        # read sst
        var2d=self.get_single_var2d(fn, 'Hwave', ifrm)
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(0,5,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    def draw2d_map_lwavep(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial peak wave period map
        '''
        varname='lwavep'
        unit='s'

        # read sst
        var2d=self.get_single_var2d(fn, 'Pwave_top', ifrm)
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(0,20,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)
        # please refer to draw2d_map_sst
        pass

    def draw2d_map_zeta(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial free surf height map
        '''
        varname='free-surf'
        unit='m'

        # read sst
        var2d=self.get_single_var2d(fn, 'zeta', ifrm)
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlWhRe
        levels=np.linspace(-3,3,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    def draw2d_map_surfcurr(self, fn, ifrm, itsk=0):
        '''
        draw 2d spatial free surf current map
        shaded contour: velocities m/s
        vectors: u,v
        '''
        varname='surfcurr'
        unit='cm/s'

        # read sst
        var2d=self.get_single_var2d(fn, 'zeta', ifrm)
        u=self.get_single_var2d(fn, 'Uwind_eastward', ifrm)
        v=self.get_single_var2d(fn, 'Vwind_northward', ifrm)
        u.values=u.values*100
        v.values=v.values*100
        var2d.values = np.power((np.power(u.values,2)+np.power(v.values,2)),0.5)
        
        tfrm=var2d['ocean_time']
        tfrm_str=tfrm.dt.strftime('%Y-%m-%d %H:%M:%S').values
        title_txt=self.roms_idom+': '+varname+' ('+unit+') '
        title_txt=title_txt+'@'+tfrm_str
        utils.write_log('TASK[%02d]: Paint %s' % (itsk, title_txt))

        ax=self.set_canvas_common()
        
        cmap=cmaps.BlGrYeOrReVi200
        levels=np.linspace(0,50,51)
        plt.contourf(
            self.lons, self.lats, 
            var2d.values,
            levels=levels, extend='both', 
            transform=ccrs.PlateCarree(), cmap=cmap)
        plt.title(title_txt, fontsize=SMFONT)

        q_mis=12 # wind vector plotting every q_miss grid
        quv = ax.quiver(self.lons[::q_mis,::q_mis],
                   self.lats[::q_mis,::q_mis],
                   u[::q_mis,::q_mis].values, v[::q_mis,::q_mis].values,
                   pivot='mid', units='inches', scale=2,
                   scale_units='inches', color="dimgray",
                   width=0.02, headwidth=3, headlength=4.5, 
                   transform=ccrs.PlateCarree())
        plt.quiverkey(quv, 1.05, -0.05, 1, r'$1 cm/s$', labelpos='N',
                       coordinates='axes')
        
        # Add a color bar
        plt.colorbar(ax=ax, shrink=0.5, extendfrac='auto')
        self.savefig(varname, tfrm.dt.strftime('%Y%m%d%H').values)

    def get_single_var2d(self, fn, varname, ifrm):
        '''get single var 2D single frame
            fn: roms file
            varname: variable name
            ifrm: frame index
        '''
        ds=xr.open_dataset(fn)
        var=ds[varname][ifrm,:,:]
        ds.close()
        return var

    def get_single_var2din3d(self, fn, varname, ifrm, ilvl=0):
        '''get single var 2D plane in 3D fields
            fn: roms file
            varname: variable name
            ifrm: frame index
            ilvl: level index
        ''' 
        ds=xr.open_dataset(fn)
        var=ds[varname][ifrm,ilvl,:,:]
        ds.close()
        return var

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
        
 
        amdn_shp=shpreader.Reader(CWD+'/shp/'+SHP_LV2).geometries()
        amdn_shp_outer=shpreader.Reader(CWD+'/shp/'+SHP_LV1).geometries()
        ax.add_geometries(
            amdn_shp_outer, ccrs.PlateCarree(),
            facecolor='none', edgecolor='black',linewidth=1, zorder = 1)
        ax.coastlines()
        
        # plot shp boundaries
        ax.add_geometries(
            amdn_shp, ccrs.PlateCarree(),
            facecolor='none', edgecolor='black',linewidth=.5, zorder = 1)
    
        return ax    
    
    def savefig(self, varname, time_str):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        fig_prefix='roms.'+self.roms_idom+'.'+varname
        
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
        plt.savefig(
                save_dir+'/'+fig_prefix+'.p%s.' % time_str, 
                dpi=120, bbox_inches='tight')
        plt.close()

