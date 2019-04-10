#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

sbatch -c10 -p largenode --mem 33G --time=7-0 --array=1-8 ./run_sl_assays.sh