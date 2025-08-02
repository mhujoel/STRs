#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

library(dplyr)
LOCUS=as.character(eval(paste(text=param[1]))); 
# ------ read in my genotypes + IBS ------#
GENO = data.table::fread(paste0( LOCUS,"_conf_het_spread.txt"))
IBS =data.table::fread(paste0("hap_neighbors.",LOCUS,".txt.gz"),h=T)
# 
IBS = IBS[IBS$ID %in% GENO$IID & IBS$IDnbr %in% GENO$IID,]
IBS_mtx = as.matrix(mutate_all(IBS %>% arrange(nbrInd) %>% 
	group_by(ID,hap) %>% 
	summarize(nbr1 = IDnbr[1], nbr1_hap = hapNbr[1], 
	nbr2 = IDnbr[2], nbr2_hap = hapNbr[2], 
	nbr3 = IDnbr[3], nbr3_hap = hapNbr[3], 
	nbr4 = IDnbr[4], nbr4_hap = hapNbr[4], 
	nbr5 = IDnbr[5], nbr5_hap = hapNbr[5],
	cM=cMlen[1], edge_cM = cMedge[1]) %>% group_by(ID) %>% mutate(n=n()) %>% 
	filter(n==2), function(x) as.numeric(as.character(x))))

# --- matrix where row is ID, column 1 is hap1 column 2 is hap2  ------ #
tab_phase = matrix(0,nrow=max(GENO$IID),ncol=2); 
tab_phase[GENO$IID,1] = GENO$A1; tab_phase[GENO$IID,2] = GENO$A2;

# iterate through individuals; update genotypes; phase using top 5 neighbors
for(ITER in 1:50){
  print(ITER)
  for (j in seq(1,nrow(IBS_mtx),by=2)){
    a1=tab_phase[IBS_mtx[j,1],1]; a2=tab_phase[IBS_mtx[j,1],2]
    
    h1=mean(c(tab_phase[IBS_mtx[j,3], IBS_mtx[j,4]],
    		tab_phase[IBS_mtx[j,5],IBS_mtx[j,6]],
    		tab_phase[IBS_mtx[j,7],IBS_mtx[j,8]],
    		tab_phase[IBS_mtx[j,9],IBS_mtx[j,10]],
    		tab_phase[IBS_mtx[j,11],IBS_mtx[j,12]]),na.rm=T)
    h2=mean(c(tab_phase[IBS_mtx[(j+1),3],IBS_mtx[(j+1),4]],
    		tab_phase[IBS_mtx[(j+1),5],IBS_mtx[(j+1),6]],
    		tab_phase[IBS_mtx[(j+1),7],IBS_mtx[(j+1),8]],
    		tab_phase[IBS_mtx[(j+1),9],IBS_mtx[(j+1),10]],
    		tab_phase[IBS_mtx[(j+1),11],IBS_mtx[(j+1),12]]),na.rm=T)
   if(h2 < h1){
      tab_phase[IBS_mtx[j,1],1]=max(a1,a2); tab_phase[IBS_mtx[j,1],2]=min(a1,a2)
     }else{
      tab_phase[IBS_mtx[j,1],1]=min(a1,a2); tab_phase[IBS_mtx[j,1],2]=max(a1,a2)
   }
  }
}
tab_phase = as.data.frame(tab_phase);  
tab_phase$IID = 1:nrow(tab_phase); 
tab_phase = tab_phase[tab_phase$IID %in% GENO$IID,]; 
colnames(tab_phase) = c("H1","H2","IID")
write.table(tab_phase[,c("IID","H1","H2")],
	paste0("phased_", LOCUS,".txt"),
	col.names=T,row.names=F,quote=F)

