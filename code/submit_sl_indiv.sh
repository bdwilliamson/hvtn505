#!/bin/bash

ml R/3.5.3-foss-2016b-fh1

sbatch -c10 --time=7-0 --array=129-448 ./run_sl_vimp.sh