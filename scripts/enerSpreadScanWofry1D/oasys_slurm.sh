#!/bin/bash -l
  
CMD_LAT=$@

CMD=$(echo /home/esrf/srio/miniconda3-py38/bin/python \
       $CMD_LAT )

echo $CMD
$CMD
