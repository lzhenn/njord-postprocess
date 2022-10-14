#!/usr/bin/env python3
'''
Postprocessing controller to dispatch all necessary
ploting scripts
This is the main script to drive the postprocessing

History:
Dec 20, 2021 --- Architecture Design, Zhenning LI

'''
import logging, logging.config
import lib, core
from utils import utils
from multiprocessing import Pool
import sys, os

CWD=sys.path[0]
def main_run():
    """Main run script"""    
    print('*************************POSTPROCESS START*************************')
    
    # logging manager
    logging.config.fileConfig(CWD+'/conf/logging_config.ini')
    
    # read log
    utils.write_log('Read Config...')
    cfg=lib.cfgparser.read_cfg(CWD+'/conf/config.fcst.post.ini')

    #------------------------------WRF Postprocessing------------------------------
    if cfg['WRF'].getboolean('wrf_flag'):
        # construnct wrf painter obj
        wrf_painter=core.wrf_painter.WRFPainter(cfg)
        # update the obj with wrf specific config paras
        wrf_painter.update(cfg)
        
        ntasks=int(cfg['WRF']['ntasks'])
         
        # --------- Plot 1D time series ---------
        if cfg['WRF'].getboolean('ts_draw'):
            wrf_painter.load_data('d03')
            stas=lib.station.construct_stas()
            stas=wrf_painter.locate_sta_pos(stas)
            wrf_painter.draw_ts_t2rh2(stas)
            wrf_painter.draw_ts_pr1h(stas)
            wrf_painter.draw_ts_wswd10(stas)
            wrf_painter.draw_ts_swslp(stas)
            wrf_painter.draw_hov_cloud(stas)
            wrf_painter.draw_hov_rhwind(stas)

        # --------- Plot 2D spatial map ---------
        
        if cfg['WRF'].getboolean('innermost_flag'):
            idom_strt=wrf_painter.wrf_num_dom
        else:
            idom_strt=1

        if cfg['WRF'].getboolean('spatial_draw'):
            for idom in range(idom_strt, wrf_painter.wrf_num_dom+1):
                dom_id='d%02d' % idom
                utils.write_log('Deal with WRFOUT %s' % dom_id)
                # load seriel output data
                wrf_painter.load_data(dom_id)
                    
                nfile=wrf_painter.wrf_num_file
                
                if ntasks > nfile:
                    ntasks=nfile
                len_per_task=nfile//ntasks

                # MULTIPROCESSING: start process pool
                process_pool = Pool(processes=ntasks)
                results=[]

                if cfg['WRF'].getboolean('debug'):
                    result=process_pool.apply_async(
                        run_mtsk, 
                        args=(0, 3, 4, wrf_painter, ))
                    results.append(result)
                    #print(results[0].get()) 
                    process_pool.close()
                    process_pool.join()    
                else:
                    # open tasks ID 0 to ntasks-1
                    for itsk in range(0, ntasks-1):  
                        
                        istrt=itsk*len_per_task    
                        iend=(itsk+1)*len_per_task-1        
                        result=process_pool.apply_async(
                            run_mtsk, 
                            args=(itsk, istrt, iend, wrf_painter, ))
                        results.append(result)

                    # open ID ntasks-1 in case of residual
                    istrt=(ntasks-1)*len_per_task
                    iend=nfile-1
                    result=process_pool.apply_async(
                        run_mtsk, 
                        args=(ntasks-1, istrt, iend, wrf_painter, ))

                    results.append(result)
                    print(results[0].get()) 
                    process_pool.close()
                    process_pool.join()    
                # MULTIPROCESSING: end process pool

                #end if: debug
            #end for: idom
        # end if: spatial draw

        # form animation
        if cfg['WRF'].getboolean('form_animation'):
            wrf_painter.form_anim(opt='d03.SST')
    
    # end if: postprocess WRF


    #------------------------------ROMS Postprocessing------------------------------
    if cfg['ROMS'].getboolean('roms_flag'):
        # construnct roms painter obj
        roms_painter=core.roms_painter.ROMSPainter(cfg)
        # update the obj with roms specific config paras
        roms_painter.update(cfg)
        ntasks=int(cfg['ROMS']['ntasks'])         
        dom_id='d02'
        utils.write_log('Deal with ROMSOUT %s' % dom_id)
        # --------- Plot 1D time series ---------
        if cfg['ROMS'].getboolean('ts_draw'):
            roms_painter.load_data(dom_id)
            stas=lib.station.construct_stas()
            stas=roms_painter.locate_sta_pos(stas)
            #roms_painter.draw_ts_sstsss(stas)
            #roms_painter.draw_ts_hwave_lwavep(stas)
            #roms_painter.draw_ts_zeta(stas)
            roms_painter.draw_hov_swt(stas, -50, 1) # sea water temperature
            roms_painter.draw_hov_sws(stas, -50, 1) # sea water salinity 

        # --------- Plot 2D spatial map ---------
        if cfg['ROMS'].getboolean('spatial_draw'):
            roms_painter.load_data(dom_id)
            nfile=roms_painter.roms_num_file
            
            if ntasks > nfile:
                ntasks=nfile
            len_per_task=nfile//ntasks

            # MULTIPROCESSING: start process pool
            process_pool = Pool(processes=ntasks)
            results=[]

            if cfg['ROMS'].getboolean('debug'):
                result=process_pool.apply_async(
                    run_mtsk_roms, 
                    args=(0, 1, 1, roms_painter, ))
                results.append(result)
                #print(results[0].get()) 
                process_pool.close()
                process_pool.join()    
            else:
                # open tasks ID 0 to ntasks-1
                for itsk in range(0, ntasks-1):  
                    
                    istrt=itsk*len_per_task    
                    iend=(itsk+1)*len_per_task-1        
                    result=process_pool.apply_async(
                        run_mtsk_roms, 
                        args=(itsk, istrt, iend, roms_painter, ))
                    results.append(result)

                # open ID ntasks-1 in case of residual
                istrt=(ntasks-1)*len_per_task
                iend=nfile-1
                result=process_pool.apply_async(
                    run_mtsk_roms, 
                    args=(ntasks-1, istrt, iend, roms_painter, ))

                results.append(result)
                #print(results[0].get()) 
                process_pool.close()
                process_pool.join()    
            # MULTIPROCESSING: end process pool

            #end if: debug
        # end if: spatial draw

        # form animation
        if cfg['ROMS'].getboolean('form_animation'):
            roms_painter.form_anim()
    # end if: postprocess ROMS

    #------------------------------SWAN Postprocessing------------------------------
    if cfg['SWAN'].getboolean('swan_flag'):
        # construnct roms painter obj
        swan_painter=core.swan_painter.SWANPainter(cfg)
        swan_painter.update(cfg)

        if cfg['SWAN'].getboolean('innermost_flag'):
            idom_strt=swan_painter.swan_num_dom
        else:
            idom_strt=1
        
        if cfg['SWAN'].getboolean('ts_draw'):
            swan_painter.load_metadata('d02')
            stas=lib.station.construct_stas()
            stas=swan_painter.locate_sta_pos(stas)
            swan_painter.draw_ts_hsig(stas)
   
        ntasks=int(cfg['SWAN']['ntasks'])
        # --------- Plot 2D spatial map ---------
        if cfg['SWAN'].getboolean('spatial_draw'):
            for idom in range(idom_strt, swan_painter.swan_num_dom+1):
                dom_id='d%02d' % idom
                        
                swan_painter.load_metadata(dom_id)
                nfile=swan_painter.swan_num_file
                
                if ntasks > nfile:
                    ntasks=nfile
                    utils.write_log('ntasks reduced to  %02d' % (ntasks))
                len_per_task=nfile//ntasks

                # MULTIPROCESSING: start process pool
                process_pool = Pool(processes=ntasks)
                results=[]

                # open tasks ID 0 to ntasks-1
                for itsk in range(0, ntasks-1):  
                    
                    istrt=itsk*len_per_task    
                    iend=(itsk+1)*len_per_task-1        
                    result=process_pool.apply_async(
                        run_mtsk_swan, 
                        args=(itsk, istrt, iend, swan_painter, ))
                    results.append(result)

                # open ID ntasks-1 in case of residual
                istrt=(ntasks-1)*len_per_task
                iend=nfile-1
                result=process_pool.apply_async(
                    run_mtsk_swan, 
                    args=(ntasks-1, istrt, iend, swan_painter, ))

                results.append(result)
                print(results[0].get()) 
                process_pool.close()
                process_pool.join()    
            # MULTIPROCESSING: end process pool
            # form animation
            if cfg['SWAN'].getboolean('form_animation'):
                swan_painter.form_anim('Hsig')

            # --------- Plot 1D time series ---------
            if cfg['SWAN'].getboolean('ts_draw'):
                swan_painter.load_data(dom_id)
   
    # end if: postprocess SWAN
    print('*************************POSTPROCESS COMPLETED*************************')

def run_mtsk_swan(itsk, istrt, iend, swan_painter):
    """
    multitask painting
    ------------------
    itask: task ID
    ilen: length of paintings per task
    painter: painter obj 
    """
    utils.write_log('TASK[%02d]: Deal with SWANOUT %s' % (
        itsk, swan_painter.swan_idom))
    # read the rest files in the list
    for ifile in range(istrt, iend+1):
        swan_painter.load_data(ifile)
        for itm in swan_painter.swan_mat_key:
            swan_painter.draw_2d_map_hsig(itm, itsk)

    return 0 



def run_mtsk_roms(itsk, istrt, iend, roms_painter):
    """
    multitask painting
    ------------------
    itask: task ID
    ilen: length of paintings per task
    painter: painter obj 
    """
    # read the rest files in the list
    for ifile in range(istrt, iend+1):
        fn=roms_painter.file_list[ifile]
        utils.write_log('TASK[%02d]: Process File: %04d of %04d --- %s' % (
                itsk, ifile, iend, fn ))
        # loop all frames in the file
        for ifrm in range(0, roms_painter.nfrms_file):
            #roms_painter.draw2d_map_sst(fn, ifrm, itsk)
            #roms_painter.draw2d_map_sss(fn, ifrm, itsk)
            roms_painter.draw2d_map_zeta(fn, ifrm, itsk)
            #roms_painter.draw2d_map_hwave(fn, ifrm, itsk)
            #roms_painter.draw2d_map_lwavep(fn, ifrm, itsk)
            #roms_painter.draw2d_map_surfcurr(fn, ifrm, itsk)
        
    return 0 


def run_mtsk(itsk, istrt, iend, wrf_painter):
    """
    multitask painting
    ------------------
    itask: task ID
    ilen: length of paintings per task
    painter: painter obj 
    """
    # read the rest files in the list
    for idx in range(istrt, iend+1):
        utils.write_log('TASK[%02d]: Process %04d of %04d --- %s' % (
                itsk, idx, iend, wrf_painter.file_list[idx]))
        wrf_painter.draw2d_map_sst(idx, itsk)
        #wrf_painter.draw2d_map_t2(idx, itsk)
        #wrf_painter.draw2d_map_rh2(idx, itsk)
        #wrf_painter.draw2d_map_wind10(idx, itsk)
        #wrf_painter.draw2d_map_slp(idx, itsk)
        #wrf_painter.draw2d_map_pr(idx, 1, itsk)
        #wrf_painter.draw2d_map_pr(idx, 3, itsk)
        #wrf_painter.draw2d_map_pr(idx, 6, itsk)
        #wrf_painter.draw2d_map_pr(idx, 12,itsk)
        #wrf_painter.draw2d_map_pr(idx, 24,itsk)
        
    return 0 


if __name__=='__main__':
    main_run()
