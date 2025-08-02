#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

library(dplyr)

tab =data.table::fread("temp.txt",h=F); 
colnames(tab) = c("IID","locus","ID","Chr","Start","End","refMotif",
	"refMotifReverseComplement","Gene","het_EUR","nAnchorIrr","nIRRpairs", "assignIRRpair",
	 "nStrictAnchorIrr","nIRRpairs_Strict","assignStrictIRRpair","nStrictAnchorIrr_H2" ,
	 "nIRRpairs_Strict_H2","assignStrictIRRpair_H2", "nStrictAnchorIrr_H3","nIRRpairs_Strict_H3",
	 "assignStrictIRRpair_H3","nStrictAnchorIrr_H4","nIRRpairs_Strict_H4","assignStrictIRRpair_H4")
	 
COVARS=read.csv("covariates.csv",h=T);
colnames(COVARS) = c("IID", "cov_SEX_GENETIC", "cov_AGE", "cov_GENETIC_ANC", paste0("ukbPC",1:10))
COVARS = COVARS[!is.na(COVARS$cov_SEX_GENETIC) & !is.na(COVARS$cov_AGE) &
	 !is.na(COVARS$cov_GENETIC_ANC) & !is.na(COVARS$ukbPC1), ]
# restrict to cov_GENETIC_ANC = 1
# Indicates samples who self-identified as 'White British' according to Field 21000 and 
# have very similar genetic ancestry based on a principal components analysis of the genotypes.
COVARS = COVARS[COVARS$cov_GENETIC_ANC == 1,]
sub_allUKB = merge(tab, COVARS,by="IID",all.y=T);

sub_allUKB$nStrictAnchorIrr_H3[is.na(sub_allUKB$nStrictAnchorIrr_H3)] = 0;
sub_allUKB$nStrictAnchorIrr_H4[is.na(sub_allUKB$nStrictAnchorIrr_H4)] = 0

# to adjust for coverage uncomment following lines and comment out last 2 commands

# COVERAGE = data.table::fread("PATH_TO_COVERAGE_DATA"); 
# sub_allUKB = merge(sub_allUKB, COVERAGE,by="IID");
# sub_allUKB$covIRR_H3 = sub_allUKB$nStrictAnchorIrr_H3/sub_allUKB$coverage; 
# sub_allUKB$covIRR_H4 = sub_allUKB$nStrictAnchorIrr_H4/sub_allUKB$coverage
# write.table(sub_allUKB[,c("IID","covIRR_H3")],"temp_h3.txt",row.names=FALSE, col.names=FALSE,quote=FALSE)
# write.table(sub_allUKB[,c("IID","covIRR_H4")],"temp_h4.txt", row.names=FALSE, col.names=FALSE,quote=FALSE)

write.table(sub_allUKB[,c("IID","nStrictAnchorIrr_H3")],"temp_h3.txt",
	row.names=FALSE, col.names=FALSE,quote=FALSE)
write.table(sub_allUKB[,c("IID","nStrictAnchorIrr_H4")],"temp_h4.txt",
	row.names=FALSE, col.names=FALSE,quote=FALSE)



