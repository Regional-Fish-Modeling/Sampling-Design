#!/bin/sh
#
#SBATCH -p normal
#SBATCH -A leetown
#SBATCH --output=%j-joboutput.out
#SBATCH -n 1
#SBATCH -c 3
#SBATCH -t 7-00:00:00
#SBATCH --array=1-2
#SBATCH --mail-user=nhitt@usgs.gov
#SBATCH --mail-type=BEGIN,END,FAIL

module load jags/4.1.0-gcc

Rscript ${SLURM_SUBMIT_DIR}/trendCovariatesSims.R

exit