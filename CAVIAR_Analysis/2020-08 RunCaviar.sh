#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=32G
#SBATCH -t 80:00:00 

module load python/3.7.4
module load caviar/1.0
module load gcc/8.3 gsl/2.3 lapack/3.7.0 blas/3.7.0
module bin caviar

python3 RunCAVIAR.py