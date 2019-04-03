#create the R script DSS.R, open nano and write the following:
setwd("/flush2/moh034/CGmap_filles")
library(DSS)
require(bsseq)
require(data.table)
OvaryT1.1 <- fread('OV-T1F1.CG_positions.min10x.tsv')
setnames(OvaryT1.1, c("chr", "pos", "N", "X"))
tail(OvaryT1.1)
dim(OvaryT1.1)
OvaryT1.2 <- fread('OV-T1F3.CG_positions.min10x.tsv')
setnames(OvaryT1.2, c("chr", "pos", "N", "X"))
OvaryT4.1 <- fread('OV-T4F3.CG_positions.min10x.tsv')
setnames(OvaryT4.1, c("chr", "pos", "N", "X"))
OvaryT4.2 <- fread('OV-T4F4.CG_positions.min10x.tsv')
setnames(OvaryT4.2, c("chr", "pos", "N", "X"))
tail(OvaryT4.2)
dim(OvaryT4.2)
BSobjovary <- makeBSseqData(list(OvaryT1.1,OvaryT1.2,OvaryT4.1,OvaryT4.2),c("OC1","OC2","OM1","OM2"))
dmlTest.sm.ov <- DMLtest(BSobjovary, group1=c("OC1","OC2"), group2=c("OM1","OM2"), smoothing =TRUE)
head(dmlTest.sm.ov)
dim(dmlTest.sm.ov)
dmls.ov <- callDML(dmlTest.sm.ov, p.threshold=0.001)
head(dmls.ov)
dim(dmls.ov)
dmrs.ov <- callDMR(dmlTest.sm.ov, delta=0.1, p.threshold=0.001)
head(dmrs.ov)
dim(dmrs.ov)
pdf("ovary-46_osbpl2dmr.pdf")
showOneDMR(dmrs.ov[46,], BSobjovary)

#create the bash script

#!/bin/bash
#SBATCH --job-name=OV-DSS_r
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mem=8gb
#SBATCH --time=6:00:00
#SBATCH --mail-user=amin.esmai@csiro.au
$SLURM_SUBMIT_DIR

module load R/3.5.0
Rscript DSS.R
