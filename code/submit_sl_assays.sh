#!/bin/bash

ml fhR/4.0.2-foss-2019b
ml jbigkit

# Args:
# 01: weighting type ("aipw" or "ipw")

sbatch -A gilbert_p -c10 --mem 33G --time=7-0 --array=1-14 ./run_sl_assays.sh ${1}
