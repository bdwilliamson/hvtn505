#!/bin/bash

ml fhR/4.0.2-foss-2019b
ml jbigkit

echo -e '#!/bin/bash \n Rscript run_sl_assays_ipw.R ' > call_run_assays.sh
chmod u+x call_run_assays.sh
sbatch -c10 --mem 33G --time=7-0 --array=1-14 ./call_run_assays.sh
rm call_run_assays.sh
