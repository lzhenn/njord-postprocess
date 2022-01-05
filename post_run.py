#!/usr/bin/env python3
'''
Postprocessing controler to dispatch all necessary
ploting scripts
This is the main script to drive the postprocessing

History:
Dec 20, 2021 --- Architecture Design, Zhenning LI

'''
import logging, logging.config
import lib, core
from utils import utils
from multiprocessing import Pool
import sys

CWD=sys.path[0]
def main_run():
    """Main run script"""    
    print('*************************POSTPROCESS START*************************')
       
    # logging manager
    logging.config.fileConfig(CWD+'/conf/logging_config.ini')
    
    # read log
    utils.write_log('Read Config...')
    cfg=lib.cfgparser.read_cfg(CWD+'/conf/config.fcst.post.ini')

    if cfg['WRF'].getboolean('wrf_flag'):
        # construnct wrf painter obj
        wrf_painter=core.wrf_painter.WRFPainter(cfg)
        # update the obj with wrf specific config paras
        wrf_painter.update(cfg)
        
        ntasks=int(cfg['WRF']['ntasks'])
        for idom in range(1, wrf_painter.wrf_num_dom+1):
            dom_id='d%02d' % idom
            utils.write_log('Deal with WRFOUT %s' % dom_id)
            # load seriel output data
            wrf_painter.load_data(dom_id)
                
            nfile=wrf_painter.wrf_num_file
            
            if ntasks > nfile:
                ntasks=nfile
            len_per_task=nfile//ntasks

            # start process pool
            process_pool = Pool(processes=ntasks)
            results=[]

            if cfg['WRF'].getboolean('debug'):
                result=process_pool.apply_async(
                    run_mtsk, 
                    args=(0, 3, 4, wrf_painter, ))
                results.append(result)
                print(results[0].get()) 
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
                #print(results[0].get()) 
                process_pool.close()
                process_pool.join()    
    print('*************************POSTPROCESS COMPLETED*************************')


def run_mtsk(itsk, istrt, iend, wrf_painter):
    """
    multitask painting
    ------------------
    itask: task ID
    ilen: length of paintings per task
    painter: painter obj 
    """
    # read the rest files in the list
    for idx in range(istrt, iend):
        utils.write_log('TASK[%02d]: Process %04d of %04d --- %s' % (
                itsk, idx, iend, wrf_painter.file_list[idx]))
        wrf_painter.draw2d_map_t2(idx, itsk)
        wrf_painter.draw2d_map_rh2(idx, itsk)
        wrf_painter.draw2d_map_wind10(idx, itsk)
        wrf_painter.draw2d_map_slp(idx, itsk)
        wrf_painter.draw2d_map_pr3h(idx, itsk)
        
    return 0 

if __name__=='__main__':
    main_run()
