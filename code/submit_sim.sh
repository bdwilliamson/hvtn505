#!/bin/bash

ml fhR/4.0.2-foss-2019b
ml jbigkit

# takes 4 command-line args:
# 01: simulation name (e.g., "aipw")
# 02: number of total Monte-Carlo reps (e.g., 1000)
# 03: the number of replicates per batch job (e.g., 5)
# 04: the i/o file prefix (e.g., "aipw")
num_n=4
njobs=`expr $2 / $3 \* $num_n`
mkdir -p $4
io_file="${4}/slurm-%A_%a.out"

echo -e "#!/bin/bash \n Rscript run_sl_assays_sim.R --sim-name ${1} --nreps-total ${2} --nreps-per-job ${3}" > ipw_sim.sh
chmod u+x ipw_sim.sh
sbatch -A gilbert_p --array=1-$njobs -e $io_file -o $io_file ./ipw_sim.sh $1 $2 $3
rm ipw_sim.sh
