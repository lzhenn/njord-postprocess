#!/usr/bin/env python3
'''
Function: load_sta_info():
Class: station 
'''
import pandas as pd
import sys
CWD=sys.path[0]

def construct_stas():
    '''construct station objs'''
    sta_df=pd.read_csv(CWD+'/db/sta.csv')
    stas=[Station(row) for row in sta_df.itertuples()]
    return stas
    
class Station():
    def __init__(self, row):
        self.id=row.Index
        self.name=row.name
        self.lat=row.lat
        self.lon=row.lon