#!/usr/bin/env python3

import os

# set SHP ROOT
SHP_ROOT= '/home/lzhenn/array74/workspace/njord_pipeline/postprocess/shp' 

try:
    os.rmdir('./fig')
    os.remove('./shp')
except:
    print('Previous link not found, fresh link')

os.system('ln -sf '+SHP_ROOT+' ./shp')
os.mkdir('fig')

