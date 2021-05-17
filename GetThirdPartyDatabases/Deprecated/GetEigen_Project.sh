#!/bin/bash
#SBATCH --account=rrg-wyeth
## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL
## CPU Usage
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-0:00
#SBATCH --nodes=1
## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error
# Load Modules
cd /project/projects/def-wyeth/DATABASES/

wget --recursive --no-parent https://xioniti01.u.hpc.mssm.edu/v1.1/

