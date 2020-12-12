#!/bin/bash

ml fhR/4.0.2-foss-2019b
ml jbigkit

echo -e '#!/bin/bash \n Rscript run_sl_assays_test_size.R --aipw ${1}' > call_run_assays_noise.sh
chmod u+x call_run_assays_noise.sh
sbatch -c10 --mem 33G --time=7-0 --array=1-50 ./call_run_assays_noise.sh ${1}
rm call_run_assays_noise.sh
