#create the R script DSS.R, open nano and write the following:
setwd("/flush2/moh034/CGmap_filles")
library(DSS)
require(bsseq)
require(data.table)
#read the ovary files and rename the headers 
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
save.image()


#read the liver files and rename the headers 

LiverT1.1 <- fread('LV-T1F3.CG_positions.min10x.tsv')
setnames(LiverT1.1, c("chr", "pos", "N", "X"))
LiverT1.2 <- fread('LV-T1F4.CG_positions.min10x.tsv')
setnames(LiverT1.2, c("chr", "pos", "N", "X"))
LiverT4.1 <- fread('LV-T4F3.CG_positions.min10x.tsv')
setnames(LiverT4.1, c("chr", "pos", "N", "X"))
LiverT4.2 <- fread('LV-T4F4.CG_positions.min10x.tsv')
setnames(LiverT4.2, c("chr", "pos", "N", "X"))

BSobjliver <- makeBSseqData(list(LiverT1.1,LiverT1.2,LiverT4.1,LiverT4.2),c("LC1","LC2","LM1","LM2"))
BSobjliver
dmlTest.sm.lv <- DMLtest(BSobjliver, group1=c("LC1","LC2"), group2=c("LM1","LM2"), smoothing =TRUE)
head(dmlTest.sm.lv)
dim(dmlTest.sm.lv)
write.table(dmlTest.sm.lv, "liver_dmlTest2.txt", sep="\t")
dmls.lv <- callDML(dmlTest.sm.lv, p.threshold=0.001)
write.table(dmlTest.sm.lv, "liver_dmls.txt", sep="\t")
head(dmls.lv)
dim(dmls.lv)
dmrs.lv <- callDMR(dmlTest.sm.lv, delta=0.1, p.threshold=0.001)
head(dmrs.lv)
dim(dmrs.lv)
write.table(dmrs.lv, "liver_dmrs2.txt", sep="\t")
pdf("liver-dmrs2.pdf")
showOneDMR(dmrs.lv[1,], BSobjliver)
dev.off() 
save.image()

#read the pituitary files and rename the headers 

PitT1.1 <- fread('Pit-T1F3.CG_positions.min10x.tsv')
setnames(PitT1.1, c("chr", "pos", "N", "X"))
PitT1.2 <- fread('Pit-T1F4.CG_positions.min10x.tsv')
setnames(PitT1.2, c("chr", "pos", "N", "X"))
PitT4.1 <- fread('Pit-T4F2.CG_positions.min10x.tsv')
setnames(PitT4.1, c("chr", "pos", "N", "X"))
PitT4.2 <- fread('Pit-T4F4.CG_positions.min10x.tsv')
setnames(PitT4.2, c("chr", "pos", "N", "X"))
BSobjpit <- makeBSseqData(list(PitT1.1,PitT1.2,PitT4.1,PitT4.2),c("PC1","PC2","PM1","PM2"))
BSobjpit
dmlTest.sm.pit <- DMLtest(BSobjpit, group1=c("PC1","PC2"), group2=c("PM1","PM2"), smoothing =TRUE)
head(dmlTest.sm.pit)
dim(dmlTest.sm.pit)
write.table(dmlTest.sm.pit, "pit_dmlTest2.txt", sep="\t")
dmls.pit <- callDML(dmlTest.sm.pit, p.threshold=0.001)
write.table(dmlTest.sm.pit, "pit_dmls.txt", sep="\t")
head(dmls.pit)
dim(dmls.pit)
dmrs.pit <- callDMR(dmlTest.sm.pit, delta=0.1, p.threshold=0.001)
head(dmrs.pit)
dim(dmrs.pit)
write.table(dmrs.pit, "pit_dmrs2.txt", sep="\t")
pdf("pit-dmrs2.pdf")
showOneDMR(dmrs.pit[1,], BSobjpit)
dev.off() 

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
