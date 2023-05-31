#!/bin/bash

while IFS=$'\t' read -r -a lines
do
  array=($lines);
  project=${array[0]};
  dir_raw=${array[1]};
  qsub -v dirPath=$dir_raw,wdir=$dir -N "beta_pval_"$project scripts/beta_pvalue.pbs;
 done < $1
