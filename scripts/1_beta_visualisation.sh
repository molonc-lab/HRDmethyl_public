#!/bin/bash

while IFS=$'\t' read -r -a lines
do
  array=($lines);
  project=${array[0]};
  qsub -v projectName=$project,wdir=$dir -N "beta_plot_"$project scripts/beta_visualisation.pbs;
 done < $1
