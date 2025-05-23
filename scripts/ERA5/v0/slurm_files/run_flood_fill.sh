#!/bin/bash

styr=1940
fnyr=2023

for yy in $( seq $styr $fnyr)
do
  sed "21s/.*/year=${yy}/" ../flood_fill.slurm > flood_fill_${yy}.slurm
  sbatch flood_fill_${yy}.slurm
done
