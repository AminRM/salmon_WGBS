#vioplot-methylation levels in bins 
install.packages('sm')
install.packages('vioplot')
library(vioplot)

setwd("C:/Users/moh034/Desktop")
require(data.table)

LC1data <- fread('LV-T1F3.CG.mbin.5mb.log', na.strings="na")
setnames(LCdata, c("chr", "pos1", "bin", "M"))
LCdatanoNA <-na.omit(LCdata)
head(LCdatanoNA)
LC1 <- LCdatanoNA$M
head(LC1)

LC2data <- fread('LV-T1F4.CG.mbin.5mb.log', na.strings="na")
setnames(LC2data, c("chr", "pos1", "bin", "M"))
LC2noNA <-na.omit(LC2data)
LC2 <- LC2noNA$M
head(LC2)

LM1data <- fread('LV-T4F3.CG.mbin.5mb.log', na.strings="na")
setnames(LM1data, c("chr", "pos1", "bin", "M"))
LM1noNA <-na.omit(LM1data)
LM1 <- LM1noNA$M
head(LM1)

LM2data <- fread('LV-T4F4.CG.mbin.5mb.log', na.strings="na")
setnames(LM2data, c("chr", "pos1", "bin", "M"))
LM2noNA <-na.omit(LM2data)
LM2 <- LM2noNA$M
head(LM2)

vioplot(LC1,LC2,LM1,LM2,names=c("LC1","LC2","LM1","LM2"), col=c("gold","gold","grey","grey"))


OC1data <- fread('OV-T1F1.CG.mbin.5mb.log', na.strings="na")
setnames(OC1data, c("chr", "pos1", "bin", "M"))
OC1datanoNA <-na.omit(OC1data)
head(OC1datanoNA)
OC1 <- OC1datanoNA$M
head(OC1)

OC2data <- fread('OV-T1F3.CG.mbin.5mb.log', na.strings="na")
setnames(OC2data, c("chr", "pos1", "bin", "M"))
OC2datanoNA <-na.omit(OC2data)
head(OC2datanoNA)
OC2 <- OC2datanoNA$M
head(OC2)


OM1data <- fread('OV-T4F3.CG.mbin.5mb.log', na.strings="na")
setnames(OM1data, c("chr", "pos1", "bin", "M"))
OM1noNA <-na.omit(OM1data)
OM1 <- OM1noNA$M
head(OM1)

OM2data <- fread('OV-T4F4.CG.mbin.5mb.log', na.strings="na")
setnames(OM2data, c("chr", "pos1", "bin", "M"))
OM2noNA <-na.omit(OM2data)
OM2 <- OM2noNA$M
head(OM2)
PC1data <- fread('Pit-T1F3.CG.mbin.5mb.log', na.strings="na")
setnames(PC1data, c("chr", "pos1", "bin", "M"))
PC1NO <-na.omit(PC1data)
PC1 <- PC1NO$M
head(PC1)
PC2data <- fread('Pit-T1F4.CG.mbin.5mb.log', na.strings="na")
setnames(PC2data, c("chr", "pos1", "bin", "M"))
PC2NO <-na.omit(PC2data)
PC2 <- PC2NO$M
head(PC2)
PM1data <- fread('Pit-T4F2.CG.mbin.5mb.log', na.strings="na")
setnames(PM1data, c("chr", "pos1", "bin", "M"))
PM1NO <-na.omit(PM1data)
PM1 <- PM1NO$M
head(PM1)
PM2data <- fread('Pit-T4F4.CG.mbin.5mb.log', na.strings="na")
setnames(PM2data, c("chr", "pos1", "bin", "M"))
PM2NO <-na.omit(PM2data)
PM2 <- PM2NO$M
head(PM2)

vioplot(OC1,OC2,OM1,OM2,names=c("OC1","OC2","OM1","OM2"), col=c("gold","gold","grey","grey"))

vioplot(PC1,PC2,PM1,PM2,names=c("PC1","PC2","PM1","PM2"), col=c("gold","gold","grey","grey"))
vioplot(PC1,PC2,LC1,LC2,names=c("PC1", "PC2","LC1","LC2"), col=c("gold","gold","grey","grey"))
vioplot(PC1,PC2,PM1,PM2,LC1,LC2,LM1,LM2, names=c("PC1", "PC2","PM1","PM2","LC1","LC2","LM1","LM2","OC1","OC2","OM1","OM2"), col=c("gold","gold","grey","grey", "gold","gold","grey","grey", "gold","gold","grey","grey"))


