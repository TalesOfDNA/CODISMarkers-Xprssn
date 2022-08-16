#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 1:00:00
#SBATCH --mail-type=FAIL,END          # Type of email notification: BEGIN,END,FAIL,A$
#SBATCH --mail-user=mayra_banuelos@brown.edu  #Email where notifications will be sent
#SBATCH -o /users/mbanuel1/data/mbanuel1/2022/CODIS_Simis/SLURM/CODIS_%A.out

module load anaconda/2020.02
source /gpfs/runtime/opt/anaconda/2020.02/etc/profile.d/conda.sh
conda activate mb_sims

#python MsPrimeSimulations_CSF1PO_CEU.py $1 $2
#python MsPrimeSimulations_D3S1358_CEU.py $1 $2
python MsPrimeSimulations_D18S51_YRI.py $1 $2


