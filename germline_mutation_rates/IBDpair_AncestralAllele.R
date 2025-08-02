#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

library(dplyr)
LOCUS=as.character(eval(paste(text=param[1]))); 
# ------ read in *phased* genotypes + IBS ------#
GENO = data.table::fread(paste0("phased_", LOCUS,".txt"),h=T)
tab_phase = matrix(0,nrow=max(GENO$IID),ncol=2); 
tab_phase[GENO$IID,1] = GENO$H1; tab_phase[GENO$IID,2] = GENO$H2;

IBS =data.table::fread(paste0("hap_neighbors.",LOCUS,".txt.gz"),h=T)
IBS = IBS[IBS$ID %in% GENO$IID & IBS$IDnbr %in% GENO$IID,]
IBS_mtx = as.matrix(mutate_all(IBS %>% arrange(nbrInd) %>% 
	group_by(ID,hap) %>% 
	summarize(nbr1 = IDnbr[1], nbr1_hap = hapNbr[1], 
	nbr2 = IDnbr[2], nbr2_hap = hapNbr[2], 
	nbr3 = IDnbr[3], nbr3_hap = hapNbr[3], 
	nbr4 = IDnbr[4], nbr4_hap = hapNbr[4], 
	nbr5 = IDnbr[5], nbr5_hap = hapNbr[5],
	nbr6 = IDnbr[6], nbr6_hap = hapNbr[6], 
	nbr7 = IDnbr[7], nbr7_hap = hapNbr[7], 
	nbr8 = IDnbr[8], nbr8_hap = hapNbr[8], 
	nbr9 = IDnbr[9], nbr9_hap = hapNbr[9], 
	nbr10 = IDnbr[10], nbr10_hap = hapNbr[10],
	cM=cMlen[1], edge_cM = cMedge[1]) %>% group_by(ID) %>% mutate(n=n()) %>% 
	filter(n==2), function(x) as.numeric(as.character(x))))

# --
phased_nbr_info = matrix(0,nrow=nrow(IBS_mtx),ncol=11)
for (j in seq(1,nrow(IBS_mtx),by=2)){
  IID = IBS_mtx[j,1]
  a1=tab_phase[IBS_mtx[j,1],1]; top_nbr_1 = tab_phase[IBS_mtx[j,3],IBS_mtx[j,4]]; 
  cM_1 =IBS_mtx[j,23]; cM_1_edge =IBS_mtx[j,24]
  a2=tab_phase[IBS_mtx[j,1],2]; top_nbr_2 = tab_phase[IBS_mtx[(j+1),3],IBS_mtx[(j+1),4]]; 
  cM_2 = IBS_mtx[(j+1),23]; cM_2_edge=IBS_mtx[(j+1),24]

  # what is ancestral allele?
  haps_1 = haps_1_id = haps_1_nbr_id =  c(); j_alt = which(IBS_mtx[,1] == IBS_mtx[j,3])[IBS_mtx[j,4]]
  for(nbr in 1:10){
  start=2*nbr+1
    haps_1 = c(haps_1 ,tab_phase[IBS_mtx[j, start],IBS_mtx[j,(start+1)]]); 
	haps_1_id = c(haps_1_id ,paste0(IBS_mtx[j, start],"_",IBS_mtx[j,(start+1)]))
	haps_1_nbr_id = c(haps_1_nbr_id ,paste0(IBS_mtx[j_alt,start],"_",IBS_mtx[j_alt,(start+1)]))
 }

  if(length(haps_1[haps_1_id %in% haps_1_nbr_id & haps_1_id != "NA_NA"]) >= 5){
    hapsCONSIDERED =  haps_1[haps_1_id %in% haps_1_nbr_id & haps_1_id != "NA_NA"]
    phased_nbr_info[j,] = c(IID,1,IBS_mtx[j,3],IBS_mtx[j,4], a1, top_nbr_1, cM_1, cM_1_edge,
    	as.numeric(names(sort(table(hapsCONSIDERED),decreasing=TRUE)[1])),
    	length(hapsCONSIDERED),as.numeric(sort(table(hapsCONSIDERED),decreasing=TRUE)[1]))
  } else{
    phased_nbr_info[j,] = c(IID,1,IBS_mtx[j,3],IBS_mtx[j,4], a1, top_nbr_1, cM_1, cM_1_edge,-9,-9,-9)
  }

  # what is ancestral allele?
  haps_2 = haps_2_id = haps_2_nbr_id  =   c(); 
  j_alt = which(IBS_mtx[,1] == IBS_mtx[(j+1),3])[IBS_mtx[(j+1),4]]
  for(nbr in 1:10){
   start=2*nbr+1
    haps_2 = c(haps_2 ,tab_phase[IBS_mtx[(j+1), start],IBS_mtx[(j+1),(start+1)]])
   haps_2_id = c(haps_2_id ,paste0(IBS_mtx[(j+1), start],"_",IBS_mtx[(j+1),(start+1)]))
    haps_2_nbr_id = c(haps_2_nbr_id , paste0(IBS_mtx[j_alt,start],"_",IBS_mtx[j_alt,(start+1)]))
 }

  if(length(haps_2[haps_2_id %in% haps_2_nbr_id  & haps_2_id != "NA_NA"]) >= 5){
    hapsCONSIDERED =  haps_2[haps_2_id %in% haps_2_nbr_id & haps_2_id != "NA_NA"]

    phased_nbr_info[(j+1),] = c(IID,2,IBS_mtx[(j+1),3],IBS_mtx[(j+1),4], a2, top_nbr_2, cM_2, cM_2_edge,
    as.numeric(names(sort(table(hapsCONSIDERED),decreasing=TRUE)[1])),length(hapsCONSIDERED),
    as.numeric(sort(table(hapsCONSIDERED),decreasing=TRUE)[1]))
  } else{
    phased_nbr_info[(j+1),] = c(IID,2, IBS_mtx[(j+1),3],IBS_mtx[(j+1),4],a2, top_nbr_2, cM_2, cM_2_edge,-9,-9,-9)
  }
}
phased_nbr_info = as.data.frame(phased_nbr_info,stringsAsFactors=F);
colnames(phased_nbr_info) = c("IID","HAP","IID_NBR","HAP_NBR","A","NBR.A","cM","cM.EDGE",
	"ANC.A","ANC.A_N_CONSIDERED","ANC.A_N")

phased_nbr_info$cM  = as.numeric(as.character(phased_nbr_info$cM)); 
phased_nbr_info$cM.EDGE  = as.numeric(as.character(phased_nbr_info$cM.EDGE))
phased_nbr_info$IID  = as.numeric(as.character(phased_nbr_info$IID)); 
phased_nbr_info$IID_NBR  = as.numeric(as.character(phased_nbr_info$IID_NBR))
phased_nbr_info$A  = as.numeric(as.character(phased_nbr_info$A)); 
phased_nbr_info$NBR.A  = as.numeric(as.character(phased_nbr_info$NBR.A))
phased_nbr_info$ANC.A  = as.numeric(as.character(phased_nbr_info$ANC.A)); 
phased_nbr_info$ANC.A_N  = as.numeric(as.character(phased_nbr_info$ANC.A_N)); 
phased_nbr_info$ANC.A_N_CONSIDERED  = as.numeric(as.character(phased_nbr_info$ANC.A_N_CONSIDERED))


phased_nbr_info_nodblcount =unique(phased_nbr_info[phased_nbr_info$cM > 5 & 
	phased_nbr_info$cM.EDGE  > 0.5 & phased_nbr_info$ANC.A_N >= 3  & 
	phased_nbr_info$ANC.A_N_CONSIDERED >= 5,]  %>% rowwise() %>% 
	mutate(ID1= ifelse(IID < IID_NBR,paste0(IID,"_",HAP),paste0(IID_NBR,"_", HAP_NBR) ) ,
		ID2 = ifelse(IID < IID_NBR,paste0(IID_NBR,"_", HAP_NBR),paste0(IID,"_",HAP))) %>%
	mutate(delta = ifelse(A == ANC.A,NBR.A- ANC.A,ifelse(NBR.A == ANC.A,A-ANC.A,-9))) %>% 
	select(ID1, ID2,cM, ANC.A,delta))

write.table(phased_nbr_info,paste0("germline_mut_", LOCUS,"_extended_info.txt"),
	col.names=T,row.names=F,quote=F)
write.table(phased_nbr_info_nodblcount,paste0("germline_mut_", LOCUS,"_unique_pairs.txt"),
	col.names=T,row.names=F,quote=F)


