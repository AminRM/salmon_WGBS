# Salmon multitissue DNA methylomes (WGBS analyses) during maturation 
Mohamed et al (in prep) Integrated transcriptome, DNA methylome and chromatin state accessibility landscapes reveal regulators of Atlantic salmon maturation

this repository contains scripts for analysing DNA methylomes from 12 WGBS libraries obtained from 3 salmon tissues (pituitary gland, ovary and liver) collected before and after initiating puberty via photoperiod manipulation. 

the following scripts were used to perform dataQC, map WGBS to the Atlantic salmon reference genome ICSASG_v2 (Lien et al., 2016), perform post-mapping QC, call methylation levels, perform differential methylation analyses using the R package DSS. 

WGBS data QC:
trimgalore_QC.sh was used to filter bases (Q scores < 30) and remove both universal and indexed adapter sequences

genome mapping:
index-bsseeker.sh was used to index the salmon genome 
map-bsseeker.sh was used to map the WGBS data to the salmon genome 

post-mapping QC: 
sortBAM.sh & pcrDUPremoval.sh were used to sort alignment files using SAMTOOLS and remove PCR duplicates using Picard tools 

methylation calling:
MethCalling.sh was used to call methylation levels (CGmap files) for all CG contexts using bs_seeker2-call_methylation.py within BSseker2 

DNA methyolme exploratory analyses:
Exploring methylomes describes exploratory analyses for the methyolme data and extraction of methylation data at CpG sites 

Differential methylation analyses:
DSS_CpG_Methylation_analysis.R for differential CpG methylation analyses using the R package DSS 

genomic distribution of DMRs:
geneScooperBed.pl and exonScooperBed.pl: perl scripts used to assign DMRs to genomic locations such as: genic (cds, intron), up-, down-stream of TSS and intergenic regions.  

Functional profiles of DMGs: 
GO_clusterProfiler.R for gene ontology (GO) enrichment among differentially methylated genes (DMGs)
