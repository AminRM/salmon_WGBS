#!/bin/bash
#SBATCH --job-name=WGBSmap
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=ALL
#SBATCH --mem=20gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR

module load python/2.7.13
module load bowtie/2.2.9/
python /flush2/moh034/BSseeker2-master/bs_seeker2-align.py -1 ${FILE}_R1.fq -2 ${FILE}_R2.fq --aligner bowtie2 --bt2-p 4 -g GCA_000233375.4_ICSASG_v2_genomic.fna -f bam -o ${FILE}.bam -u ${FILE}_unmapped
