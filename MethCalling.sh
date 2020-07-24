
#!/bin/bash
#SBATCH --job-name=MCall
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR

module load python/2.7.13
python /flush2/moh034/BSseeker2-master/bs_seeker2-call_methylation.py -i ${FILE}_filtered.bam --sorted -o ${FILE} --db /flush2/moh034/BSseeker2-master/bs_utils/reference_genomes/GCA_000233375.4_ICSASG_v2_genomic.fna_bowtie2


