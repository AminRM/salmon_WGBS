#!/bin/bash
#SBATCH --job-name=OV-T1F1-MS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mem=12gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR
module load samtools/1.9.0
samtools sort -@ ${SLURM_NTASKS} -T temp -O BAM ${FILE}.bam -o ${FILE}_sorted.bam
