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
        
        # load seriel output data
        wrf_painter.load_data('d01')

        # paint frame 0 t2m
        wrf_painter.draw2d_map_t2(0)
    
    print('*************************POSTPROCESS COMPLETED*************************')

if __name__=='__main__':
    main_run()