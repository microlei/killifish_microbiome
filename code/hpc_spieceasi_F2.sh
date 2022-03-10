#!/bin/bash
#SBATCH --partition=bigmem           # Queue selection
#SBATCH --job-name=spieceasi        # Job name
#SBATCH --mail-type=ALL               # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lma@whoi.edu    # Where to send mail
#SBATCH --ntasks=4                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task (if 8 == 8 threads)
#SBATCH --mem=1000000                     # Job memory request (if 10,000 == 10000MB == 10GB)
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=hpc_logs/spieceasi_F2.log   # Standard output/error
#SBATCH --qos=unlim
#export OMP_NUM_THREADS=10       # 8 threads

#conda activate r4
source activate spieceasi
pwd; hostname; date

Rscript hpc_spieceasi_F2.R

date
