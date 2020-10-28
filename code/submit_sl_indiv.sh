#!/bin/bash

ml fhR/4.0.2-foss-2019b
ml jbigkit

## arg 1 is the risk type

sbatch -c10 --mem 33G --time=7-0 --array=191-610 ./run_sl_vimp.sh $1
