#!/bin/bash

# Report the total node memory every 60 seconds

./check_mem.sh -t 60 > checkmem-${SLURMD_NODENAME}-${SLURM_JOB_ID}.out
