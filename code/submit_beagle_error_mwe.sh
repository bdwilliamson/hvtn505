#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

sbatch ./beagle_error_mwe.sh

sbatch -M beagle ./beagle_error_mwe.sh