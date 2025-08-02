#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

LOCUS=as.character(eval(paste(text=param[1]))); REP=as.character(eval(paste(text=param[2])));

COVARS=read.csv("covariates.csv",h=T);
colnames(COVARS) = c("IID", "cov_SEX_GENETIC", "cov_AGE", "cov_GENETIC_ANC", paste0("ukbPC",1:10))
COVARS = COVARS[!is.na(COVARS$cov_SEX_GENETIC) & !is.na(COVARS$cov_AGE) &
	 !is.na(COVARS$cov_GENETIC_ANC) & !is.na(COVARS$ukbPC1), ]
# restrict to cov_GENETIC_ANC = 1
# Indicates samples who self-identified as 'White British' according to Field 21000 and 
# have very similar genetic ancestry based on a principal components analysis of the genotypes.
COVARS = COVARS[COVARS$cov_GENETIC_ANC == 1,]

chr2=data.table::fread("chr2.raw",h=T); chr2=chr2[,c("IID","2:47416318G:A_G")]
colnames(chr2) = c("IID","MSH2")
chr5=data.table::fread("chr5.raw",h=T); chr5=chr5[,c("IID","5:80638411G:T_G")]
colnames(chr5) = c("IID","MSH3")

tab =data.table::fread("temp.txt",h=F); 
colnames(tab) = c("IID","locus","ID","Chr","Start","End","refMotif",
	"refMotifReverseComplement","Gene","het_EUR","nAnchorIrr","nIRRpairs", "assignIRRpair", 
	"nStrictAnchorIrr","nIRRpairs_Strict","assignStrictIRRpair","nStrictAnchorIrr_H2" ,
	"nIRRpairs_Strict_H2","assignStrictIRRpair_H2", "nStrictAnchorIrr_H3","nIRRpairs_Strict_H3",
	"assignStrictIRRpair_H3","nStrictAnchorIrr_H4","nIRRpairs_Strict_H4","assignStrictIRRpair_H4")

imputed_h3 = data.table::fread(paste0("IRRs_phased_imp.",LOCUS,".",REP,".h3.txt.gz"));
colnames(imputed_h3) = c("ID","IRRs_H3","hap1phased_H3","hap2phased_H3","hap1imp_H3","hap2imp_H3"); 
imputed_h3$imputedH3 = imputed_h3$hap1imp_H3 + imputed_h3$hap2imp_H3
imputed_h4 = data.table::fread(paste0("IRRs_phased_imp.",LOCUS,".",REP,".h4.txt.gz")); 
colnames(imputed_h4) = c("ID","IRRs_H4","hap1phased_H4","hap2phased_H4","hap1imp_H4","hap2imp_H4"); 
imputed_h4$imputedH4 = imputed_h4$hap1imp_H4 + imputed_h4$hap2imp_H4

sub_allUKB = merge(tab, COVARS,by="IID",all.y=T); 
sub_allUKB = merge(sub_allUKB, chr2,by="IID"); sub_allUKB = merge(sub_allUKB, chr5,by="IID"); 
sub_allUKB = merge(sub_allUKB, imputed_h3,by.x="IID",by.y="ID");  
sub_allUKB = merge(sub_allUKB, imputed_h4,by.x="IID",by.y="ID");

# to adjust for coverage uncomment following lines (and end of command to compute H3 & H4)
# COVERAGE = data.table::fread("PATH_TO_COVERAGE_DATA"); 
# sub_allUKB = merge(sub_allUKB, COVERAGE,by="IID");
# sub_allUKB$covIRR_H3 = sub_allUKB$nStrictAnchorIrr_H3/sub_allUKB$coverage; 
# sub_allUKB$covIRR_H4 = sub_allUKB$nStrictAnchorIrr_H4/sub_allUKB$coverage

sub_allUKB$nStrictAnchorIrr_H3[is.na(sub_allUKB$nStrictAnchorIrr_H3)] = 0;  
sub_allUKB$nStrictAnchorIrr_H4[is.na(sub_allUKB$nStrictAnchorIrr_H4)] = 0
sub_allUKB$H3 = sub_allUKB$nStrictAnchorIrr_H3 # /sub_allUKB$coverage; 
sub_allUKB$H4 = sub_allUKB$nStrictAnchorIrr_H4 # /sub_allUKB$coverage; 

mdl_ukb_h3 = lm(sub_allUKB,formula = H3~imputedH3); 
sub_allUKB$residH3 = sub_allUKB$H3 - predict(mdl_ukb_h3,sub_allUKB);  
cropVAL_h3 = sd(sub_allUKB$residH3)* sqrt(nrow(sub_allUKB)/100); 
sub_allUKB$residH3_cropped = pmax(pmin(sub_allUKB$residH3, cropVAL_h3),-cropVAL_h3)

mdl_ukb_h4 = lm(sub_allUKB,formula = H4~imputedH4);
sub_allUKB$residH4 = sub_allUKB$H4 - predict(mdl_ukb_h4,sub_allUKB); 
cropVAL_h4 = sd(sub_allUKB$residH4)* sqrt(nrow(sub_allUKB)/100); 
sub_allUKB$residH4_cropped = pmax(pmin(sub_allUKB$residH4, cropVAL_h4),-cropVAL_h4)

regression_h3 = summary(lm(sub_allUKB,
	formula=residH3_cropped~MSH2+MSH3+cov_AGE+ 
	ukbPC1+ukbPC2+ukbPC3+ukbPC4+ukbPC5+ukbPC6+ukbPC7+ukbPC8+ukbPC9+ukbPC10 + 
	cov_SEX_GENETIC))$coefficient[2:4,]

regression_h4 = summary(lm(sub_allUKB,
	formula=residH4_cropped~MSH2+MSH3+cov_AGE+ 
	ukbPC1+ukbPC2+ukbPC3+ukbPC4+ukbPC5+ukbPC6+ukbPC7+ukbPC8+ukbPC9+ukbPC10 + 
	cov_SEX_GENETIC))$coefficient[2:4,]

ALL = data.frame(LOCUS = LOCUS, REP = REP, DEF = c("H3","H4"), 
	Nnonzero = c(sum(sub_allUKB$nStrictAnchorIrr_H3 > 0), 
	sum(sub_allUKB$nStrictAnchorIrr_H4 > 0)) ,
	MSH2beta = c(regression_h3[1,1], regression_h4[1,1]) , 
	MSH2se = c(regression_h3[1,2], regression_h4[1,2]) , 
	MSH2p = c(regression_h3[1,4], regression_h4[1,4]) ,
	 MSH3beta = c(regression_h3[2,1], regression_h4[2,1]) , 
	 MSH3se = c(regression_h3[2,2], regression_h4[2,2]) , 
	 MSH3p = c(regression_h3[2,4], regression_h4[2,4]) , 
	 AGEbeta = c(regression_h3[3,1], regression_h4[3,1]) , 
	 AGEse = c(regression_h3[3,2], regression_h4[3,2]) , 
	 AGEp = c(regression_h3[3,4], regression_h4[3,4]))

# add some additional columns that will be useful for GWAS:
sub_allUKB$cov_AGE_SQ = sub_allUKB$cov_AGE^2
sub_allUKB$imputedH3_cov_AGE = sub_allUKB$imputedH3 * sub_allUKB$cov_AGE
sub_allUKB$imputedH4_cov_AGE = sub_allUKB$imputedH4 * sub_allUKB$cov_AGE
sub_allUKB$imputedH3_cov_AGE_SQ = sub_allUKB$imputedH3 * sub_allUKB$cov_AGE_SQ
sub_allUKB$imputedH4_cov_AGE_SQ = sub_allUKB$imputedH4 * sub_allUKB$cov_AGE_SQ
sub_allUKB$FID = sub_allUKB$IID

sub_allUKB = sub_allUKB[,c("FID","IID", "locus", "ID", "Chr", "Start", "End", "refMotif", 
	"refMotifReverseComplement","Gene", "het_EUR", "nAnchorIrr", "nIRRpairs", "assignIRRpair",
	"nStrictAnchorIrr", "nIRRpairs_Strict", "assignStrictIRRpair","nStrictAnchorIrr_H2",
	 "nIRRpairs_Strict_H2", "assignStrictIRRpair_H2","nStrictAnchorIrr_H3", "nIRRpairs_Strict_H3",
	  "assignStrictIRRpair_H3","nStrictAnchorIrr_H4", "nIRRpairs_Strict_H4", "assignStrictIRRpair_H4",
	  "cov_SEX_GENETIC", "cov_AGE", "cov_GENETIC_ANC", "ukbPC1", "ukbPC2","ukbPC3", "ukbPC4", 
	  "ukbPC5", "ukbPC6", "ukbPC7", "ukbPC8", "ukbPC9","ukbPC10", "MSH2", "MSH3", "IRRs_H3",
	   "hap1phased_H3", "hap2phased_H3","hap1imp_H3", "hap2imp_H3", "imputedH3", "IRRs_H4",
	    "hap1phased_H4","hap2phased_H4", "hap1imp_H4", "hap2imp_H4", "imputedH4", "H3",
	    "H4", "residH3", "residH3_cropped", "residH4", "residH4_cropped","cov_AGE_SQ",
	     "imputedH3_cov_AGE", "imputedH4_cov_AGE", "imputedH3_cov_AGE_SQ","imputedH4_cov_AGE_SQ")]
	     
write.table(sub_allUKB,paste0("phenotypes.",LOCUS,".",REP,".UKB.txt"),
	row.names=FALSE, col.names=T,quote=FALSE)

write.table(ALL,paste0("associations.",LOCUS,".",REP,".UKB.txt"),
	row.names=FALSE, col.names=T,quote=FALSE)
