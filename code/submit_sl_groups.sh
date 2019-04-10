#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

sbatch -c4 --time=7-0 --array=1-128 ./run_sl_vimp.sh