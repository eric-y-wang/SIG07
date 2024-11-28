#!/bin/bash

# The name of the job
#BSUB -J glm
# The name of the queue you wan't to submit the job to
#BSUB -q cpuqueue
# The number of processors you want
#BSUB -n 32
# Sets memory requirements for the job (per processor)
#BSUB -R rusage[mem=20]
# How long the job will take (you job will be killed if it runs over the time limit)
#BSUB -W 48:00
# Output and error log files (optional but recommended)
#BSUB -o /lila/data/rudensky/EYW/R-out.%J
#BSUB -e /lila/data/rudensky/EYW/R-err.%J

# load bulkseq conda environment
source ~/.bashrc
mamba activate R-deseq

# set directory (with fail safe in case it fails)
cd /lila/data/rudensky/EYW/git_projects/SIG07/analysis|| { echo "Failure"; exit 1; }

R interaction_glmGamPoi.r
