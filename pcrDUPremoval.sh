#!/bin/bash
#SBATCH --job-name=DUPrem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load picard
java -Xmx8g -jar /flush2/moh034/picard.jar MarkDuplicates I=${FILE}_sorted.bam O=${FILE}_filtered.bam M=${FILE}_dup_metrics.txt REMOVE_DUPLICATES=true AS=true SORTING_COLLECTION_SIZE_RATIO=0.10
