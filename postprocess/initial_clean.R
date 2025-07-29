#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

suppressWarnings(suppressMessages(library(dplyr))); library(stringr); 
options(dplyr.summarise.inform = FALSE); library(stringdist)
ID=as.character(eval(paste(text=param[1])));

# read in relevant tables:
strict_defn = data.table::fread("motif_2to6bp_correspondence_perfectIRR_noRepeat.txt.gz",h=T);
polymorphicSTR = data.table::fread(paste0("polyEUR_2to6bp_noHomopolymer_STRs_",ID,".txt"),h=T)
reads = data.table::fread(paste0("reads_sub_",ID,".txt"),h=T)

# how many repeats do we require to consider a sequence potentially of that motif:
polymorphicSTR = merge(polymorphicSTR, 
	data.frame(repCt=c(1,2,3,4,5,6),expected= c(145,70,45,32,25,20)),by="repCt")
	
polymorphicSTR$locus = paste0(polymorphicSTR$Chrom,":",round(polymorphicSTR$Start/1e+6,digits=1))
reads = reads %>% rowwise() %>% mutate(locus = paste0(RNEXT,":",round(PNEXT/1e+6,digits=1)) )
# ------- 
gt_tab = merge(reads[abs(reads$POS - reads$PNEXT) >= 2 | reads$RNEXT != reads$RNAME,c("locus", "QNAME", "SEQ","RNEXT","PNEXT")],
	polymorphicSTR[,c("Motif", "ID", "Chrom", "Start", "End", "Gene","het_EUR","locus","refMotif","refMotifReverseComplement","expected")],
	by=c("locus")) %>%
	rowwise() %>% mutate(delta = abs(PNEXT - Start)) %>% filter(delta < 5e+3) %>% 
	rowwise() %>% 
	mutate(matches = max(c(stringr::str_count(SEQ,refMotif),stringr::str_count(SEQ,refMotifReverseComplement))) - expected, 
		which = c(refMotif,refMotifReverseComplement)[which.max(c(stringr::str_count(SEQ,refMotif),stringr::str_count(SEQ,refMotifReverseComplement)))] ) %>% 
		filter(matches > 0) %>% group_by(QNAME) %>% 
		mutate(select = which.min(delta),id =seq(1,n())) %>% 
		filter(select==id) %>% rowwise() %>% 
		mutate( strictIRR = ifelse( sum(stringr::str_detect(SEQ, strict_defn$perfectrep[strict_defn$relevantMotif == which])) > 0,1,0) ,
		  hamming_dist = min(stringdist(SEQ, strict_defn$perfectrep_long[strict_defn$relevantMotif == which], method="hamming")) ) %>%  
		  mutate(strictIRR_H2 = ifelse(hamming_dist <= 2,1,0),strictIRR_H3 = ifelse(hamming_dist <= 3,1,0) ,strictIRR_H4 = ifelse(hamming_dist <= 4,1,0)  ) %>% 
		  group_by(locus,ID,Chrom,Start,End,refMotif,refMotifReverseComplement,Gene,het_EUR) %>%
		   mutate(mates=paste0(unique(RNEXT),":",min(PNEXT),"-",max(PNEXT)))
gt_tab$IID = ID

# --- how to determine motif for IRR pairs:
reads_dup = data.table::fread(paste0("reads_pair_sub_",ID,".txt"),h=T)
assignment_df = merge( data.frame(unique_rep = unique(c(gt_tab$refMotif, gt_tab$refMotifReverseComplement)), reps = stringr::str_length(unique(c(gt_tab$refMotif, gt_tab$refMotifReverseComplement))) ), data.frame(reps=c(1,2,3,4,5,6),expected= c(145,70,45,32,25,20))); assignment_df$unique_rep = as.character(assignment_df$unique_rep)
rotations = unique(data.frame(Motif = c(gt_tab$refMotif, gt_tab$refMotifReverseComplement), refMotif =  c(gt_tab$refMotif, gt_tab$refMotif), refMotifReverseComplement =  c(gt_tab$refMotifReverseComplement, gt_tab$refMotifReverseComplement)))

duplicate_info = reads_dup %>% rowwise() %>% mutate(Motif = as.character(assignment_df$unique_rep[which.max(stringr::str_count(SEQ, assignment_df$unique_rep) - assignment_df$expected)])) %>% rowwise() %>% mutate(matches_seen =  stringr::str_count(SEQ,Motif) - assignment_df$expected[assignment_df$unique_rep==Motif] ) %>% filter(matches_seen > 0) %>% rowwise() %>% mutate(strictIRR = ifelse( sum(stringr::str_detect(SEQ, strict_defn$perfectrep[strict_defn$relevantMotif == Motif])) > 0,1,0) , 
hamming_dist = min(stringdist(SEQ, strict_defn$perfectrep_long[strict_defn$relevantMotif == Motif],method="hamming")) ) %>% inner_join(rotations,by="Motif") %>% group_by(QNAME,refMotif,refMotifReverseComplement) %>% summarize(nReads = n(),nReads_strict = sum(strictIRR), nReads_strict_H2 = sum(hamming_dist <= 2), nReads_strict_H3 = sum(hamming_dist <= 3), nReads_strict_H4 = sum(hamming_dist <= 4) ) %>% filter(nReads == 2) %>% group_by(refMotif,refMotifReverseComplement) %>% summarize(nIRRpairs = n(),nIRRpairs_Strict = sum(nReads_strict==2) ,nIRRpairs_Strict_H2 = sum(nReads_strict_H2==2),nIRRpairs_Strict_H3 = sum(nReads_strict_H3==2),nIRRpairs_Strict_H4 = sum(nReads_strict_H4 ==2))

# ---------
gt_tab = merge(gt_tab, duplicate_info,by=c("refMotif","refMotifReverseComplement"),all.x=T); gt_tab$nIRRpairs[is.na(gt_tab$nIRRpairs)] = 0; gt_tab$nIRRpairs_Strict[is.na(gt_tab$nIRRpairs_Strict)] = 0;  gt_tab$nIRRpairs_Strict_H2[is.na(gt_tab$nIRRpairs_Strict_H2)] = 0; gt_tab$nIRRpairs_Strict_H3[is.na(gt_tab$nIRRpairs_Strict_H3)] = 0; gt_tab$nIRRpairs_Strict_H4[is.na(gt_tab$nIRRpairs_Strict_H4)] = 0;

# ---------
write.table(gt_tab[,c( "IID","locus","mates","QNAME","RNEXT","PNEXT", "strictIRR","strictIRR_H2","strictIRR_H3","strictIRR_H4","ID", "Chrom", "Start", "End", "Gene","refMotif","refMotifReverseComplement",  "Motif","het_EUR","nIRRpairs","nIRRpairs_Strict","nIRRpairs_Strict_H2","nIRRpairs_Strict_H3","nIRRpairs_Strict_H4")], paste0("initial_clean_",ID,".txt"),row.names=F,col.names=T,quote=F,sep="\t");
