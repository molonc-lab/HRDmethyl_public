#!/bin/bash

qsub -v working_dir=${dir} ${dir}/scripts/snp_Cosmic_signature.pbs

