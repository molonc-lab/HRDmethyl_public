#!/bin/bash

while IFS=$'\t' read -r -a lines
do
  array=($lines);
  project=${array[0]};
  dir_raw=${array[1]};
  
  #only methylation detected tumour type had the mRNA file downloaded
  mRNA_file=$(ls ${dir_raw}/mRNA | wc -l)
  if (( $mRNA_file > 0 ));
    then
      qsub -v projectName=$project,dirPath=$dir_raw,wdir=$dir -N "mRNA_"$project scripts/mRNA_edgeR.pbs;
  fi
 done < $1
