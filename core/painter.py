#!/usr/bin/env python3
'''
Base Class: Painter
'''

import sys

from utils import utils

print_prefix='lib.painter>>'
CWD=sys.path[0]

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
        