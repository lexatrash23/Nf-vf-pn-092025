#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # Partition to submit to
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o venom_test_%j.out    # File to which STDOUT will be written, including job ID

# need to load the java & python modules
module load jdk
module load python

# optionally, load your nextflow conda environment using the conda
# source activate nextflow

nextflow run main.nf -with-conda -profile cannon \
    -with-timeline results/timeline.html \
    -with-report results/report.html \
    -with-trace results/trace.txt \
    -with-dag results/flowchart.png
