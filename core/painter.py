#!/usr/bin/env python3
'''
Base Class: Painter
'''

import sys

from utils import utils
import matplotlib.pyplot as plt
import os

print_prefix='lib.painter>>'
CWD=sys.path[0]

# Constants
BIGFONT=22
MIDFONT=18
SMFONT=14

class Painter:

    '''
    Construct core class painter
    
    Attributes
    -----------
    
    Methods
    -----------
        __init__

    '''
    def __init__(self, cfg):
        utils.write_log(print_prefix+'construct painter')
        self.arch_root=cfg['INPUT']['arch_root']
        self.fig_root=cfg['OUTPUT']['fig_root']
        self.init_ts=cfg['INPUT']['model_init_ts']
 
    def set_y1y2_axis(self,sta):
        '''set common xy plot style, 
            including left and right y axis and common x axis'''

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

    def savefig_ts(self, figname):
        '''save fig portal'''
        save_dir=self.fig_root+'/'+self.init_ts
        if not(os.path.isdir(save_dir)):
                os.mkdir(save_dir)        
        plt.savefig(
                save_dir+'/'+figname, 
                dpi=120, bbox_inches='tight')
        plt.close()

