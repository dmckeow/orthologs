#!/bin/bash
#SBATCH --no-requeue
#SBATCH --mem 16G
#SBATCH -p genoa64
#SBATCH --qos normal
#SBATCH -t 12:00:00
#SBATCH --output=logs/slurm-misc.%j.out
#SBATCH --error=logs/slurm-misc.%j.err

to_delete="/users/asebe/dmckeown/projects/crg-bcaortho/work/"
find $to_delete -type f -delete

find $to_delete -type d -delete
