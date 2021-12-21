#/usr/bin/env python
"""Commonly used utilities
    Function    
    ---------------
    throw_error(msg):
    write_log(msg, lvl=20):
"""
import sys
import logging

CWD=sys.path[0]

def throw_error(msg):
    '''
    throw error and exit
    '''
    logging.error(msg)
    exit()

def write_log(msg, lvl=20):
    '''
    write logging log to log file
    level code:
        CRITICAL    50
        ERROR   40
        WARNING 30
        INFO    20
        DEBUG   10
        NOTSET  0
    '''
    logging.log(lvl, msg)