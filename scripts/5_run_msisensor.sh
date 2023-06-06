#!/bin/bash
while read -a line
do
  echo ${line[0]}
  qsub -v normal=${line[3]},tumour=${line[2]},sampleid=${line[0]} -N MSI_SENSOR_${line[0]} run_msisensor.pbs  
done < $1
