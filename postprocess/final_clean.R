#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)
ID=as.character(eval(paste(text=param[1])));

suppressWarnings(suppressMessages(library(dplyr))); options(dplyr.summarise.inform = FALSE)
init = data.table::fread(paste0("initial_clean_",ID,".txt"),h=T); 
mates = data.table::fread(paste0(ID,"_mates_quality.txt"),h=F)

# assign IRRs + mapping quality of mates

combo  = merge(init, mates, by.x=c("QNAME","RNEXT","PNEXT"), by.y=c("V1","V2","V3"))  %>% 
	filter(V4 >= 30) %>% 
	group_by(IID,locus,ID,Chrom,Start,End,refMotif,refMotifReverseComplement,Gene,het_EUR,nIRRpairs,nIRRpairs_Strict,nIRRpairs_Strict_H2,nIRRpairs_Strict_H3,nIRRpairs_Strict_H4) %>% 
	summarize(nAnchorIrr=n(),nStrictAnchorIrr = sum(strictIRR),
		nStrictAnchorIrr_H2= sum(strictIRR_H2),nStrictAnchorIrr_H3 = sum(strictIRR_H3),
		nStrictAnchorIrr_H4 = sum(strictIRR_H4)) %>% filter(nAnchorIrr > 0) %>% 
		group_by(refMotif) %>% 
		mutate(assignIRRpair = ifelse(n() > 1 & mean(nIRRpairs) > 0,0,1), 
		assignStrictIRRpair = ifelse(sum(nStrictAnchorIrr > 0) > 1  & mean(nIRRpairs_Strict) > 0,0,1) ,
		 assignStrictIRRpair_H2 = ifelse(sum(nStrictAnchorIrr_H2 > 0) > 1  & mean(nIRRpairs_Strict_H2) > 0,0,1) ,
		  assignStrictIRRpair_H3 = ifelse(sum(nStrictAnchorIrr_H3 > 0) > 1  & mean(nIRRpairs_Strict_H3) > 0,0,1) ,
		   assignStrictIRRpair_H4 = ifelse(sum(nStrictAnchorIrr_H4 > 0) > 1  & mean(nIRRpairs_Strict_H4) > 0,0,1) )

combo$Gene = stringr::str_remove_all(combo$Gene," ")

write.table(combo[,c("IID", "locus", "ID", "Chrom", "Start", "End", "refMotif","refMotifReverseComplement", "Gene", "het_EUR", "nAnchorIrr","nIRRpairs", "assignIRRpair", "nStrictAnchorIrr","nIRRpairs_Strict","assignStrictIRRpair", "nStrictAnchorIrr_H2","nIRRpairs_Strict_H2","assignStrictIRRpair_H2", "nStrictAnchorIrr_H3","nIRRpairs_Strict_H3","assignStrictIRRpair_H3", "nStrictAnchorIrr_H4","nIRRpairs_Strict_H4","assignStrictIRRpair_H4")],paste0("clean_",ID,".txt"),
	row.names=F,col.names=T,quote=F);
