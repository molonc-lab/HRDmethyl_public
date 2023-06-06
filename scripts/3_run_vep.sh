#!/bin/bash

qsub -v working_dir=${dir} ${dir}/scripts/vep_customCV_gnomADe.pbs

