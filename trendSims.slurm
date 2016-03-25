#!/bin/sh


#SBATCH -p normal
#SBATCH -A leetown
#SBATCH --output=%j-joboutput.out
#SBATCH -n 1
#SBATCH -c 3
#SBATCH -t 7-00:00:00
#SBATCH --array=1-2
#SBATCH --mail-user=echildress@usgs.gov
#SBATCH --mail-type=END,FAIL

module load jags/4.1.0-intel

Rscript ${SLURM_SUBMIT_DIR}/trendSims.R --args $SLURM_ARRAYID %j

exit
