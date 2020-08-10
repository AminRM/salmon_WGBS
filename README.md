# Salmon DNA methylomes (WGBS analyses) during maturation 
Mohamed et al (2020) Integrated transcriptome, DNA methylome and chromatin state accessibility landscapes reveal regulators of Atlantic salmon maturation

this repository contains scripts for WGBS data QC, genome mapping, methylation calling using BSseeker2, ethylome exploratory analyses and differential CpG methylation from WGBS data using the R package DSS. 

the following scripts were used to perform dataQC, map WGBS to the Atlantic salmon reference genome ICSASG_v2 (Lien et al., 2016), perform post-mapping QC, call methylation levels, perform differential methylation analyses. 

1- trimgalore_QC.sh was used to filter bases (Q scores < 30) and remove both universal and indexed adapter sequences

2- index-bsseeker.sh was used to index the salmon genome 

3- map-bsseeker.sh was used to map the WGBS data to the salmon genome 

4- sortBAM.sh & pcrDUPremoval.sh were used to sort alignment files using SAMTOOLS and remove PCR duplicates using Picard tools 

5- MethCalling.sh was used to call methylation levels (CGmap files) for all CG contexts using bs_seeker2-call_methylation.py within BSseker2 

6- Exploring methylomes describes exploratory analyses for the methyolme data and extraction of methylation data at CpG sites 

7- DSS_CpG_Methylation_analysis.R for differential CpG methylation analyses using the R package DSS 

8- geneScooperBed.pl and exonScooperBed.pl: perl scripts used to assign DMRs to genomic locations such as: genic (cds, intron), up-, down-stream of TSS and intergenic regions.  

9- GO_clusterProfiler.R for gene ontology (GO) enrichment among differentially methylated genes (DMGs)
