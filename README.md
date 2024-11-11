nextflow run main.nf -profile local

sbatch submit_nf.sh main.nf -profile slurm