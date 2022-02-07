#!/usr/bin/env python3
'''
Batch run controller for Postprocessing 
to dispatch all necessary ploting scripts
            Zhenning LI
History:
Jan 27, 2022 --- Build
'''
import lib
import sys, os

CWD=sys.path[0]
def main_run():

    cfg=lib.cfgparser.read_cfg(CWD+'/conf/config.fcst.post.ini.smp')
    dir_info=os.listdir(cfg['INPUT']['arch_root'])
    for info in dir_info:
        cfg['INPUT']['model_init_ts']=info
        print('++++++++++++++++++++++++++BATCH RUN: %s++++++++++++++++++++++++++' % info)
        lib.cfgparser.write_cfg(cfg, CWD+'/conf/config.fcst.post.ini')
        os.system('python3 '+CWD+'/post_run.py')

if __name__=='__main__':
    main_run()
