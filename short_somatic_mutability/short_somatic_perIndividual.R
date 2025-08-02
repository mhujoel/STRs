#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

# Get command line arguments
ID=as.character(eval(paste(text=param[1]))); 
DAT=as.character(eval(paste(text=param[2]))); 

# Load necessary libraries without printing startup messages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

# Extract start and end repeat sequences, and flank positions from DAT
STARTstr = stringr::str_split(DAT,"_")[[1]][5]
repLen=stringr::str_length(stringr::str_split(DAT,"_")[[1]][7])
STARTstr = stringr::str_split(DAT,"_")[[1]][5]; 
ENDstr = stringr::str_split(DAT,"_")[[1]][6]
L_FLANK = as.numeric(stringr::str_split(stringr::str_remove(stringr::str_split(DAT,"_")[[1]][8],"L"),",")[[1]])
R_FLANK = as.numeric(stringr::str_split(stringr::str_remove(stringr::str_split(DAT,"_")[[1]][9],"R"),",")[[1]])

# Load individual read data
tab=data.table::fread(paste0("IID_",ID,".txt"),h=F,sep=" "); 

# Builds a consensus sequence from multiple reads and associated qualities
# Aligns sequences by a known start location and calculates the most frequent base at each position
consensus = function(SEQUENCES,QUALITIES,START){
  if(length(SEQUENCES) == 1){ return(SEQUENCES) }
  else{
    START_LOCS = as.data.frame(stringr::str_locate(string= SEQUENCES, START)) ; SEQUENCES = stringr::str_split(SEQUENCES,"");
    QUALITIES = stringr::str_split(QUALITIES,""); SEQdf = data.frame(SPOT=seq(-156,156,by=1),SEQ="",QUAL="",stringsAsFactors = F)
    for(i in 1:length(SEQUENCES)){
      # first REP location start is at START_LOCS[i,2]-2 ; center REP to start at "0"
      ADD = data.frame(SPOT=seq(1,length(SEQUENCES[[i]]))-START_LOCS[i,2]-2, SEQ=SEQUENCES[[i]],QUAL=QUALITIES[[i]])
      SEQdf[SEQdf$SPOT %in% ADD$SPOT,"SEQ"] = paste0(SEQdf[SEQdf$SPOT %in% ADD$SPOT,"SEQ"],ADD$SEQ);
      SEQdf[SEQdf$SPOT %in% ADD$SPOT,"QUAL"] = paste0(SEQdf[SEQdf$SPOT %in% ADD$SPOT,"QUAL"],ADD$QUAL)
    }
    SEQdf = SEQdf[SEQdf$SEQ != "",]; CONSENSUS = NULL
    for(i in 1:nrow(SEQdf)){
      # restrict to high quality bases:
      opts= stringr::str_split(SEQdf$SEQ[i],"")[[1]][utf8ToInt(SEQdf$SEQ[i])-33 >= 25]
      # restrict to high quality bases:
      # opts =  stringr::str_split(SEQdf$SEQ[i],"")[[1]][stringr::str_split(SEQdf$QUAL[i],"")[[1]] %in% c("F",":") ]
      if(length(opts) > 0){
        sortedVAL = sort(table(opts),decreasing=TRUE)
        if(length(sortedVAL)==1){
          CONSENSUS = paste0(CONSENSUS,names(sortedVAL[1]),collapse = "")
        }else{
          if(sortedVAL[1] == sortedVAL[2]){ # if a tie --> store X
            CONSENSUS = paste0(CONSENSUS,"X",collapse="")
          }else{
            CONSENSUS = paste0(CONSENSUS,names(sortedVAL[1]),collapse = "")
          }
        }
      }else{
        CONSENSUS = paste0(CONSENSUS,"X",collapse="")
      }
    }
  }
  return(CONSENSUS)
}

mismatch_bases <- function(strand,jump,som_seq,som_qual,seq,start,end,L_flank_set,R_flank_set,repeatLength){
  # Check base quality at key bases used to match flanks. 
  locs_flanks = as.data.frame(stringr::str_locate(som_seq,c(start,end)));
  if(nrow(locs_flanks) != 2 | sum(is.na(locs_flanks)) > 0){
    # didnt have flanks!
    L_FLANK_GOOD = NA; R_FLANK_GOOD = NA; L_FLANK_SUB_GOOD =NA; R_FLANK_SUB_GOOD=NA
    return(data.table::data.table(rateSTRT_mis = NA, nSTRT_mis = NA,
                                  flnk_l_hiP = NA,flnk_l_hiP_sub = NA,
                                  flnk_r_hiP = NA, flnk_r_hiP_sub = NA,  
                                  LEN =NA,Q = NA,NONF_SCORE =NA,
                                  HI_QUAL_Q =NA,HI_QUAL_LEN =NA,NONF_HIQUAL = NA))
  }else{
    L_FLANK = stringr::str_split(som_qual,"")[[1]][locs_flanks[1,1]:locs_flanks[1,2]]
    R_FLANK = stringr::str_split(som_qual,"")[[1]][locs_flanks[2,1]:locs_flanks[2,2]]
    L_FLANK_SUB = L_FLANK[L_flank_set[!is.na(L_flank_set)]]; R_FLANK_SUB = R_FLANK[R_flank_set[!is.na(R_flank_set)]]
    
    L_FLANK_GOOD = mean(unlist(lapply(L_FLANK,FUN = function(x) utf8ToInt(x)-33)) >= 25)
    R_FLANK_GOOD = mean(unlist(lapply(R_FLANK,FUN = function(x) utf8ToInt(x)-33)) >= 25)
    L_FLANK_SUB_GOOD = mean(unlist(lapply(L_FLANK_SUB,FUN = function(x) utf8ToInt(x)-33)) >= 25)
    R_FLANK_SUB_GOOD = mean(unlist(lapply(R_FLANK_SUB,FUN = function(x) utf8ToInt(x)-33)) >= 25)
    
    if(L_FLANK_GOOD == "NaN"){L_FLANK_GOOD = NA}; if(L_FLANK_SUB_GOOD == "NaN"){L_FLANK_SUB_GOOD = NA}
    if(R_FLANK_GOOD == "NaN"){R_FLANK_GOOD = NA}; if(R_FLANK_SUB_GOOD == "NaN"){R_FLANK_SUB_GOOD = NA}
    
  }
  
  # SEQ mismatch score for the flank at the beginning of the read. 
  # Currently we check SEQ and QUAL at the end of the read to look for PCR error. 
  # It now strikes me that misalignments will often cause SEQ to mismatch at the beginning of the read,
  # so I think we could easily filter most such errors by computing (1) score = mismatch rate; and 
  # (2) # of base pairs checked. For this check I'd restrict to reasonable-quality bases (QUAL>=':'),
  # for which I think we'd expect to see near-0 mismatches.
  if(strand == "reverse"){ 
    locs_E = as.data.frame(stringr::str_locate(seq,end)); 
    locs_E_som = as.data.frame(stringr::str_locate(som_seq,end));
    if(is.na(locs_E[1,1])){
      START_MISMATCH_RATE=START_MISMATCH_BP_CHECKED=-9
    }else{
      ref = stringr::str_split(seq,"")[[1]][locs_E[1,2]:stringr::str_length(seq)]
      compar = stringr::str_split(som_seq,"")[[1]][locs_E_som[1,2]:stringr::str_length(som_seq)]
      compar_qual = stringr::str_split(som_qual,"")[[1]][locs_E_som[1,2]:stringr::str_length(som_seq)]
      set = pmin(length(ref),length(compar))
      compar=compar[1:set];ref= ref[1:set]; 
      compar_qual=compar_qual[1:set]
      compar_qual_values = unlist(lapply(compar_qual,FUN = function(x) utf8ToInt(x)-33))
     START_MISMATCH_RATE = mean(compar[compar_qual_values>=25 & ref != "X"] != ref[compar_qual_values>=25  & ref != "X"])
      START_MISMATCH_BP_CHECKED = sum(compar_qual_values>=25 & ref != "X")
     }
  }else{
    locs_S = as.data.frame(stringr::str_locate(seq,start)); 
    if(is.na(locs_S[1,1])){
      START_MISMATCH_RATE=START_MISMATCH_BP_CHECKED=-9
    }else{
      locs_S_som = as.data.frame(stringr::str_locate(som_seq,start));
      ref = rev(stringr::str_split(seq,"")[[1]][1:locs_S[1,1]])
      compar = rev(stringr::str_split(som_seq,"")[[1]][1:locs_S_som[1,1]])
      compar_qual = rev(stringr::str_split(som_qual,"")[[1]][1:locs_S_som[1,1]])
      set = pmin(length(ref),length(compar))
      compar=compar[1:set];ref= ref[1:set]; compar_qual=compar_qual[1:set]
      compar_qual_values = unlist(lapply(compar_qual,FUN = function(x) utf8ToInt(x)-33))
      START_MISMATCH_RATE = mean(compar[compar_qual_values>=25 & ref != "X"] != ref[compar_qual_values>=25  & ref != "X"])
      START_MISMATCH_BP_CHECKED = sum(compar_qual_values>=25 & ref != "X")
    }
  }
  
  if(strand == "reverse"){ 
    # for reverse reads --> find start flank; align there and predict mismatch based on reference sequence
    locs = as.data.frame(stringr::str_locate(seq,start));
    if(is.na(locs[1,1])){# didn't have perfect end flank! (e.g. may be an X)
      sub_som_qual = NA;sub_som_qual_hi=NA
    } else{ 
      sub_seq_none = stringr::str_split(seq,"")[[1]][1:(locs[1,2]-repeatLength)]; # flank with no REP 
      sub_seq_one = stringr::str_split(seq,"")[[1]][1:(locs[1,2])]; # flank with one REP
      sub_seq_two = stringr::str_split(seq,"")[[1]][1:(locs[1,2]+repeatLength)]; # flank with two REP
      # from the view point of the somatic read: 
      locs_som = as.data.frame(stringr::str_locate(som_seq,start));
      if(jump < 0){ # if jump - --> want to start at "no REP"
        if(jump == -1){
          pred_mismatch = which(rev(sub_seq_none)  != rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) &
                                  rev(sub_seq_none) != "X" & rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) != "X")
          pred_match = which(rev(sub_seq_none)  == rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) &
                               rev(sub_seq_none) != "X" & rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) != "X")
        }
        if(jump == -2){
          pred_mismatch = which(rev(sub_seq_none)  != rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) & 
                                  rev(sub_seq_none) != "X" &  rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) != "X")
          pred_match = which(rev(sub_seq_none)  == rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) & 
                               rev(sub_seq_none) != "X" &  rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) != "X")
        }
        sub_som_qual = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2]-repeatLength)])[pred_mismatch]; # flank with no "GCA"
        sub_som_qual_hi = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2]-repeatLength)])[pred_match]; # flank with no "GCA"
      } 
      # if jump + --> want to consider the "jumped" REP as potential mismatch
      if(jump == 1){
        pred_mismatch = which(rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) != rev(sub_seq_none) &
                                rev(sub_seq_none) != "X" & rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) != "X")
        pred_match = which(rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) == rev(sub_seq_none) &
                             rev(sub_seq_none) != "X" & rev(sub_seq_one[(repeatLength+1):length(sub_seq_one)]) != "X")
        sub_som_qual = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2])])[pred_mismatch]; # flank with one "GCA"
        sub_som_qual_hi = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2])])[pred_match]; # flank with one "GCA"
      }
      if(jump == 2){
        pred_mismatch = which(rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) != rev(sub_seq_none) & 
                                rev(sub_seq_none) != "X" & rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)])  != "X")
        pred_match = which(rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)]) == rev(sub_seq_none) & 
                             rev(sub_seq_none) != "X" & rev(sub_seq_two[(2*repeatLength+1):length(sub_seq_two)])  != "X")
        sub_som_qual = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2]+repeatLength)])[pred_mismatch]; # flank with two "GCA"
        sub_som_qual_hi = rev(stringr::str_split(som_qual,"")[[1]][1:(locs_som[1,2]+repeatLength)])[pred_match]; # flank with two "GCA"
      }
    }}
  if(strand == "forward"){ 
    # for forward reads --> find end flank; align there and predict mismatch based on reference sequence
    locs = as.data.frame(stringr::str_locate(seq,end));end_len = stringr::str_length(seq)
    if(is.na(locs[1,1])){# didn't have perfect end flank! (e.g. may be an X)
      sub_som_qual = NA;sub_som_qual_hi=NA
    } else{
      sub_seq_none = stringr::str_split(seq,"")[[1]][(locs[1,1]+repeatLength):end_len]; # flank with no REP
      sub_seq_one = stringr::str_split(seq,"")[[1]][(locs[1,1]):end_len]; # flank with one REP
      sub_seq_two = stringr::str_split(seq,"")[[1]][(locs[1,1]-repeatLength):end_len]; # flank with two REP
      # from the view point of the somatic read: 
      locs_som = as.data.frame(stringr::str_locate(som_seq,end));som_end=stringr::str_length(som_qual) 
      if(jump < 0){ # if jump - --> want to start at "no REP"
        if(jump == -1){
          pred_mismatch = which(sub_seq_none != sub_seq_one[1:length(sub_seq_none)] & 
                                  sub_seq_none != "X" &  sub_seq_one[1:length(sub_seq_none)] != "X")
          pred_match = which(sub_seq_none == sub_seq_one[1:length(sub_seq_none)] & 
                               sub_seq_none != "X" &  sub_seq_one[1:length(sub_seq_none)] != "X")
          
        }
        if(jump == -2){
          pred_mismatch = which(sub_seq_none != sub_seq_two[1:length(sub_seq_none)] & 
                                  sub_seq_none != "X" & sub_seq_two[1:length(sub_seq_none)] != "X")
          pred_match = which(sub_seq_none == sub_seq_two[1:length(sub_seq_none)] & 
                               sub_seq_none != "X" & sub_seq_two[1:length(sub_seq_none)] != "X")
        }
        
        sub_som_qual = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]+repeatLength):som_end][pred_mismatch]; # flank with no "GCA"
        sub_som_qual_hi = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]+repeatLength):som_end][pred_match]; # flank with no "GCA"
      } 
      # if jump + --> want to consider the "jumped" REP as potential mismatch
      if(jump == 1){
        pred_mismatch = which(sub_seq_one[1:length(sub_seq_none)] != sub_seq_none & 
                                sub_seq_none != "X" & sub_seq_one[1:length(sub_seq_none)] != "X")
        pred_match = which(sub_seq_one[1:length(sub_seq_none)] == sub_seq_none & 
                             sub_seq_none != "X" & sub_seq_one[1:length(sub_seq_none)] != "X")
        sub_som_qual = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]):som_end][pred_mismatch]; # flank with one "GCA"
        sub_som_qual_hi = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]):som_end][pred_match]; # flank with one "GCA"
      }
      if(jump == 2){
        pred_mismatch = which(sub_seq_two[1:length(sub_seq_none)] != sub_seq_none & 
                                sub_seq_none != "X" & sub_seq_two[1:length(sub_seq_none)] != "X")
        pred_match = which(sub_seq_two[1:length(sub_seq_none)] == sub_seq_none & 
                             sub_seq_none != "X" & sub_seq_two[1:length(sub_seq_none)] != "X")
        sub_som_qual = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]+repeatLength):som_end][pred_mismatch]; # flank with two "GCA"
        sub_som_qual_hi = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]+repeatLength):som_end][pred_match]; # flank with two "GCA"
      }
    }
  }
  # qualities at potential mismatch bases: paste0(sub_som_qual[!is.na(sub_som_qual)],collapse="")
  # mismatch at potential mismatch bases: baseDISAGREEMENT = seq_1[1:subset] != seq_2[1:subset] 
  # qualities at predicted hi-quality (no mismatch bases): paste0(sub_som_qual_hi[!is.na(sub_som_qual_hi)],collapse="")
  
  sub_som_qual_values = unlist(lapply(sub_som_qual[!is.na(sub_som_qual)],FUN = function(x) utf8ToInt(x)-33)) 
  sub_som_qual_hi_values = unlist(lapply(sub_som_qual_hi[!is.na(sub_som_qual_hi)],FUN = function(x) utf8ToInt(x)-33)) 
  
  # starting sequence:  START_MISMATCH_RATE ; START_MISMATCH_BP_CHECKED
  # starting flank: L_FLANK_GOOD, R_FLANK_GOOD, L_FLANK_SUB_GOOD, R_FLANK_SUB_GOOD
  return(data.table::data.table(rateSTRT_mis = START_MISMATCH_RATE, nSTRT_mis = START_MISMATCH_BP_CHECKED,
                                flnk_l_hiP = L_FLANK_GOOD,flnk_l_hiP_sub = L_FLANK_SUB_GOOD,
                                flnk_r_hiP = R_FLANK_GOOD, flnk_r_hiP_sub = R_FLANK_SUB_GOOD,  
                                LEN = length(sub_som_qual[!is.na(sub_som_qual)]),
                                Q = paste0(sub_som_qual[!is.na(sub_som_qual)],collapse=""),
                                NONF_SCORE = mean(sub_som_qual_values < 30),
                                HI_QUAL_Q = paste0(sub_som_qual_hi[!is.na(sub_som_qual_hi)],collapse=""),
                                HI_QUAL_LEN = length(sub_som_qual_hi[!is.na(sub_som_qual_hi)]),
                                NONF_HIQUAL = mean(sub_som_qual_hi_values < 30) 
  ))
}

# Group reads by individual and allele, generate consensus, count reads
# Only keep individuals with at least two distinct repeat lengths, each supported by 3+ reads
# Check difference in length between alleles is meaningful (>= 5 repeat units)
tab_sum_temp = tab %>% group_by(V1,V3) %>% 
  mutate(nREADs=n(),
         midSEG = names(sort(table(as.character(V4)),decreasing=T)[1]),
         midSEGlen = stringr::str_length(names(sort(table(as.character(V4)),decreasing=T)[1])),
         SEQ = consensus(SEQUENCES= as.character(V5),QUALITIES= as.character(V6),START= STARTstr),
         QUAL = ifelse(n() == 1, as.character(V6),"Consensus")) %>% 
  group_by(V1) %>% mutate(nALLELES=length(unique(midSEGlen))) %>% filter(nALLELES >= 2) %>% group_by(V1, nALLELES) %>% 
  mutate( A1len = midSEGlen[order(nREADs,decreasing=T)][1], 
          A2len = midSEGlen[order(nREADs,decreasing=T)][midSEGlen[order(nREADs,decreasing=T)] != A1len][1],
          A1consensus = SEQ[order(nREADs,decreasing=T)][1], 
          A2consensus = SEQ[order(nREADs,decreasing=T)][midSEGlen[order(nREADs,decreasing=T)] != A1len][1],
          A1midconsensus = midSEG[order(nREADs,decreasing=T)][1], 
          A2midconsensus = midSEG[order(nREADs,decreasing=T)][midSEGlen[order(nREADs,decreasing=T)] != A1len][1]) %>% 
  group_by(V1) %>%
  filter(nREADs[order(nREADs,decreasing=T)][midSEGlen[order(nREADs,decreasing=T)] != A1len][1] >= 3 & 
           sum(nREADs[match(unique(midSEGlen), midSEGlen)] > 2) == 2 & 
           abs(unique(A1len) -unique(A2len))  >= 5*repLen) 
if(nrow(tab_sum_temp) > 0){
  # Identify reads that don't exactly match either of the two main alleles but are within ±2 repeat units
  # Compare against consensus sequences to calculate "jump"
  # For each jump type (-2, -1, 1, 2), call mismatch_bases()
  # Use quality and mismatch filters to define somatic events
tab_sum = tab_sum_temp %>% ungroup() %>% rowwise() %>% 
  mutate(SOM_POT = ifelse(abs(midSEGlen-A1len) <= 2*repLen | abs(midSEGlen-A2len) <= 2*repLen,1,0)) %>% filter(SOM_POT==1) %>%
  rowwise() %>% 
  mutate(fromLEN = ifelse(midSEGlen==A1len | midSEGlen==A2len, NA, c(A1len,A2len)[which.min(c(abs(midSEGlen-A1len),abs(midSEGlen-A2len)))]),
         fromSEQ = ifelse(midSEGlen==A1len , A1consensus,ifelse(midSEGlen==A2len, A2consensus,c(as.character(A1consensus),as.character(A2consensus))[which.min(c(abs(midSEGlen-A1len), abs(midSEGlen-A2len)))]))) %>%
  mutate(jump = ifelse(is.na(fromLEN),NA,(midSEGlen-fromLEN)/repLen)) %>% filter(jump %in% c(-2,-1,1,2,NA)) %>%
  rowwise() %>%
  mutate(j_neg2 =  mismatch_bases(strand=as.character(V2),jump=-2,
                                  som_seq=as.character(V5),som_qual=as.character(V6),
                                  seq=as.character(fromSEQ),start= STARTstr,end=ENDstr,
                                  L_flank_set=L_FLANK,R_flank_set=R_FLANK,repeatLength=repLen),
         j_neg1 =  mismatch_bases(strand=as.character(V2),jump=-1,
                                  som_seq=as.character(V5),som_qual=as.character(V6),
                                  seq=as.character(fromSEQ),start= STARTstr,end=ENDstr,
                                  L_flank_set=L_FLANK,R_flank_set=R_FLANK,repeatLength=repLen),
         j_pos1 =  mismatch_bases(strand=as.character(V2),jump=1,
                                  som_seq=as.character(V5),som_qual=as.character(V6),
                                  seq=as.character(fromSEQ),start= STARTstr,end=ENDstr,
                                  L_flank_set=L_FLANK,R_flank_set=R_FLANK,repeatLength=repLen),
         j_pos2 =  mismatch_bases(strand=as.character(V2),jump=2,
                                  som_seq=as.character(V5),som_qual=as.character(V6),
                                  seq=as.character(fromSEQ),start= STARTstr,end=ENDstr,
                                  L_flank_set=L_FLANK,R_flank_set=R_FLANK,repeatLength=repLen)) %>%
  rowwise() %>% mutate(
    somatic_neg2 = ifelse(j_neg2$LEN >= 4 & j_neg2$NONF_SCORE < 0.2 & j_neg2$rateSTRT_mis < 0.05 & 
                             ((j_neg2$flnk_l_hiP_sub == 1 | is.na(j_neg2$flnk_l_hiP_sub)) & 
                                (j_neg2$flnk_r_hiP_sub == 1 | is.na(j_neg2$flnk_r_hiP_sub))),1,0),
    somatic_neg1 = ifelse(j_neg1$LEN >= 4 & j_neg1$NONF_SCORE < 0.2 & j_neg1$rateSTRT_mis < 0.05 & 
                            ((j_neg1$flnk_l_hiP_sub == 1 | is.na(j_neg1$flnk_l_hiP_sub)) & 
                               (j_neg1$flnk_r_hiP_sub == 1 | is.na(j_neg1$flnk_r_hiP_sub))),1,0),
    somatic_pos1 = ifelse(j_pos1$LEN >= 4 & j_pos1$NONF_SCORE < 0.2 & j_pos1$rateSTRT_mis < 0.05 & 
                            ((j_pos1$flnk_l_hiP_sub == 1 | is.na(j_pos1$flnk_l_hiP_sub)) & 
                               (j_pos1$flnk_r_hiP_sub == 1 | is.na(j_pos1$flnk_r_hiP_sub))),1,0),
    somatic_pos2 = ifelse(j_pos2$LEN >= 4 & j_pos2$NONF_SCORE < 0.2 & j_pos2$rateSTRT_mis < 0.05 & 
                            ((j_pos2$flnk_l_hiP_sub == 1 | is.na(j_pos2$flnk_l_hiP_sub)) & 
                               (j_pos2$flnk_r_hiP_sub == 1 | is.na(j_pos2$flnk_r_hiP_sub))),1,0)) %>%
   select(c(V1,V3,V4,jump,fromLEN,A1len,A1midconsensus,A2len,A2midconsensus,
            somatic_neg2,somatic_neg1,somatic_pos1,somatic_pos2))
# --- save data
write.table(tab_sum,paste0("summary_somatic_",ID,".txt"),row.names=F,col.names=T,quote=F)
} 
