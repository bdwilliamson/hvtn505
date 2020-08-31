#!/bin/bash

ml R/3.5.3-foss-2018b
ml jbigkit

sbatch -c10 --mem 33G --time=7-0 --array=1-14 ./run_sl_assays.sh
