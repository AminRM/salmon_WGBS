#!/bin/bash
#SBATCH --job-name=2index-ref
#SBATCH --output  dump.out
#SBATCH --error dump.err
#SBATCH --mail-type=ALL
#SBATCH --mem=20gb
#SBATCH --time=16:00:00
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR

module load python
module load bowtie/2.2.9/
python BSseeker2-master/bs_seeker2-build.py -f GCA_000233375.4_ICSASG_v2_genomic.fna --aligner bowtie2
