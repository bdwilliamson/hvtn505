#!/bin/bash

ml fhR/4.0.2-foss-2019b

sbatch -c10 --mem 33G --time=7-0 ./run_sl.sh
