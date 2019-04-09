#!/bin/bash

ml R/3.5.3-foss-2016b-fh1

sbatch -c10 -p largenode --mem 33G --time=7-0 --array=1-128 ./run_sl_vimp.sh