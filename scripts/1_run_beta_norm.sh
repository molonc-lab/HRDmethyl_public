#!/bin/bash

while IFS=$'\t' read -r -a lines
do
  array=($lines);
  project=${array[0]};
  dir_raw=${array[1]};
  qsub -v projectName=$project,dirPath=$dir_raw,wdir=$dir -N "norm_beta_"$project scripts/beta_norm.pbs;
 done < $1
