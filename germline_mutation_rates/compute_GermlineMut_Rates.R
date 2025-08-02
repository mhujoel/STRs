#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

library(dplyr)

LOCUS=as.character(eval(paste(text=param[1]))); print(LOCUS)
ADJUSTMENT_FACTOR=as.numeric(eval(paste(text=param[2])));

tab = read.table(paste0("germline_mut_", LOCUS,"_unique_pairs.txt"),h=T)
# -------- compute TMRCA
Ne = read.table("autosome10k_ibdne.ne.txt",h=T)
cM_to_TMRCA  <- function(cM,Ne_vec){
    IBDlen = cM / 100; # convert to Morgans
    genMax = 500; # numerically integrate out to genMax generations
    numer = 0; denom = 0; pNotCoal = 1;
    for (gen in 1:genMax){
        pCoal = 1 / (2*Ne_vec[min(gen,length(Ne_vec))]); # assume constant Ne beyond Ne_vec
        thisSegCoal = pCoal * pNotCoal * (4*gen^2) * exp(-2*gen*IBDlen);
        numer = numer + gen*thisSegCoal; denom = denom + thisSegCoal;
        pNotCoal = pNotCoal * (1-pCoal);
    }
    exp_TMRCA = numer / denom; return(exp_TMRCA)
}

tab$TMRCA = cM_to_TMRCA(tab$cM , Ne_vec= Ne$NE[-1]); tab$gen = 2 * tab$TMRCA

# -------- compute rates for common alleles
df_return = tab %>% filter(abs(delta) <= 2) %>% group_by(ANC.A) %>%
	 mutate(totalGEN = sum(gen),totalN = n()) %>% group_by(ANC.A,totalGEN,totalN,delta) %>% 
	 summarize(n=n()) %>% group_by(ANC.A) %>% mutate(rate = n/totalGEN) %>% 
	 filter(totalN >= 500)
 # ----- add back missing jumps
toADD = unique(df_return[,c("ANC.A","totalGEN","totalN")])
df_return = df_return %>% filter(delta != 0)
tomerge = rbind( cbind(toADD,delta=-2),cbind(toADD,delta=-1),
				cbind(toADD,delta=1),cbind(toADD,delta=2))
df_return = merge(df_return,tomerge,all=T)
df_return$n[is.na(df_return$n)] = 0;  df_return$rate[is.na(df_return$rate)] = 0;
df_return$L_exact = df_return$U_exact = -9
for(i in 1:nrow(df_return)){
sum = binom.test(x=df_return$n[i],n=round(df_return$totalGEN[i]))$conf.int; 
df_return$L_exact[i] =sum[1]; df_return$U_exact[i] = sum[2]
}
df_return$rate = df_return$rate * ADJUSTMENT_FACTOR
df_return$L_exact = df_return$L_exact * ADJUSTMENT_FACTOR
df_return$U_exact = df_return$U_exact * ADJUSTMENT_FACTOR
# ---------
write.table(df_return,paste0("germline_mut_rates_",LOCUS,".txt"),
	col.names=T,row.names=F,quote=F)