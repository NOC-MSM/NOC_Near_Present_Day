#!/bin/bash

styr=1940
fnyr=2023

for yy in $( seq $styr $fnyr)
do
  sed "21s/.*/year=${yy}/" ../flood_fill_waves_tier2.slurm > flood_fill_waves_tier2_${yy}.slurm
  sbatch flood_fill_waves_tier2_${yy}.slurm
done
