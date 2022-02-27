#!/bin/bash

# call mkmod in parallel

# NOTE: split reg_5x5 file first
# ----
#    split -l 100 --numeric-suffixes=100 -a 3 reg_5x5 reg_5x5_
#    split -l 150 -d -a 1 reg_5x5 reg_5x5_
#    ./mkmod.sh year

# define year
year=$1

# number of processors
NP=15

for ((i=101; i<101+$NP; ++i))
do
  #nohup Rscript mkmod.r $year $i $NP > log/nohup.out$i &
  echo "./mkmod reg_5x5_$i $year > nohup.out$i &"
  nohup ./mkmod reg_5x5_$i $year > nohup.out$i &
  sleep 3s
done

echo "Making surface data in the background mode..."
