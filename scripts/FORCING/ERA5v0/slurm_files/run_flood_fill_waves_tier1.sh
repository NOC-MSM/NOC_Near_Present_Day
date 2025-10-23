#!/bin/bash

styr=1940
fnyr=2023

for yy in $( seq $styr $fnyr)
do
  sed "21s/.*/year=${yy}/" ../flood_fill_waves_tier1.slurm > flood_fill_waves_tier1_${yy}.slurm
  sbatch flood_fill_waves_tier1_${yy}.slurm
done
