#!/bin/bash

# call mkmod in parallel

# NOTE: split reg_5x5 file first
# ----
#    split -l 150 --numeric-suffixes=100 -a 3 reg_5x5 reg_5x5_
#    split -l 150 -d -a 1 reg_5x5 reg_5x5_
#    ./mkmod.sh

# number of processors
NP=10

for ((i=0; i<$NP; ++i))
do
  #nohup Rscript mkmod.r 2005 $i $NP > log/nohup.out$i &
  echo "./mkmod reg_5x5_$i 2005 > nohup.out$i &"
  nohup ./mkmod reg_5x5_$i 2005 > nohup.out$i &
  sleep 3s
done

echo "Making surface data in the background mode..."
