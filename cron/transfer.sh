#!/bin/sh
# This script is routinly run to fetch the operational
# forecast result from the COAWST system.

SRC=/home/pathop/njord/data/archive_calypso/
AIM=/home/metctm1/array/data/pathop/
INIT_DIR=`date -d "-2 days" "+%Y%m%d"`
mkdir -p $AIM/${INIT_DIR}00
scp -P 18425 pathop@124.88.34.133:${SRC}/${INIT_DIR}00/hsign_* ${AIM}/${INIT_DIR}00/
scp -P 18425 pathop@124.88.34.133:${SRC}/${INIT_DIR}00/hswell_* ${AIM}/${INIT_DIR}00/

# Clean previous data and home ncx files
# clean outdated archived data
CLEAN_DATE=`date -d "-14 days" "+%Y%m%d"`
CLEAN_DIR=${AIM}/${CLEAN_DATE}00
rm -rf $CLEAN_DIR