#!/bin/bash -l
#
# Usage ./slurm_submit.sh CONFIGURATION_FILE(S)
#
if [ $# -eq 0 ]
then
   echo "Error: no filename given"
   exit 1
fi
# https://confluence.esrf.fr/display/SCKB/Partitions
FILE=$@
echo $FILE
#sbatch -p nice --nodes=1 --mincpus=28 --time=12:00:00 --job-name=$FILE \
#sbatch -p nice-long --nodes=1 --mincpus=28 --time=99:59:00 --job-name=$FILE \
sbatch -p nice --nodes=1 --mincpus=28 --time=12:00:00 --job-name=$FILE \
   /home/esrf/srio/paper-undulators-resources/scripts/enerSpreadScanWofry1D/oasys_slurm.sh \
   $FILE
