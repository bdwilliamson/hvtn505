#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

sbatch -M beagle -c10 -p largenode --mem 33G --time=7-0 --array=129-448 ./run_sl_vimp.sh