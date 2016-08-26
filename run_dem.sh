#!/bin/bash

cd /tmp

CHUNK_SIZE=512
ARRAY_NAME="aia2_${CHUNK_SIZE}"
DEM_NAME="dem_${CHUNK_SIZE}"

  
echo -e "Start"
date
/usr/bin/time -f "Elapsed Time: %E" iquery -aqf "store(dem(${ARRAY_NAME}),${DEM_NAME})"
echo -e "Stop"
date
echo -e "-------------------"
more dem_log*
