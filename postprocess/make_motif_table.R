#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

READ_LEN=as.numeric(eval(paste(text=param[1])));

# 2-6 BP; want all possible "circular" versions + reverse complement 
suppressWarnings(suppressMessages(library(dplyr))); options(dplyr.summarise.inform = FALSE)
install.packages("gtools");
options_seq = NULL
# create all potential 2-6 bp repeat motifs
for(i in 2:6){
  options_seq = c(options_seq,
                  apply(gtools::permutations(n=4,r=i,v=c("A","C","T","G"),
                                             repeats.allowed=TRUE),1,function(x) paste0(x,collapse=""))
  )
}
options_seq = unique(options_seq)

# string equivalence function; see whether 2 strings are "circular" versions of themselves
# e.g. AGC and CAG
ff = function(x, y) (nchar(y) == nchar(x)) && (grepl(y, strrep(x, 2), fixed = TRUE))

# populate data frame of potential motifs:
all = data.frame(Motif = as.character(options_seq[1]), 
                 refMotif = as.character(options_seq[1]), 
                 revComp = 0)

for(i in 2:length(options_seq)){
  where = as.character(all$refMotif); s = as.character(options_seq[i]); 
  potentialMatch  = reverseCompliment = NULL
  
  for(j in 1:length(where)){ 
    # see if a circular version of sequence (or reverse complement) is already in list as reference motifs
    circular_version = ff(where[j],s)
    circular_version_reversecomplement = ff(where[j],stringi::stri_reverse(chartr("ATGC","TACG",s)))
    # if a potential match has occurred, note it + add if its reverse complement
    potentialMatch = c(potentialMatch, circular_version | circular_version_reversecomplement )
    reverseCompliment = c(reverseCompliment, circular_version_reversecomplement )
  }
  
  if(mean(potentialMatch==FALSE) < 1 ){
    # there is a matching reference motif 
    if(mean(reverseCompliment==FALSE) < 1){
      # the sequence is the reverse complement of a "known" reference motif
      all = rbind(all, data.frame(Motif = s, refMotif = where[potentialMatch][1], revComp = 1))
    } else{
      all = rbind(all, data.frame(Motif = s, refMotif = where[potentialMatch][1], revComp = 0))
    }
  } else{
    # never seen this motif before --> add it as  new reference motif
    all = rbind(all, data.frame(Motif = s, refMotif = s,  revComp = 0))
  }
}
# 
all$refMotif = as.character(all$refMotif); all$Motif = as.character(all$Motif); 

# fix repeat motifs of length 4 and 6 that could be duplicates (or triplicates) 
# of repeat motifs of length 2 and 3; e.g. CGCG = CG
correctRef = data.frame(refMotif_true =  unique(all$refMotif[stringr::str_length(all$refMotif) %in% c(2,3)])) %>% 
  rowwise() %>% mutate(refMotif_dbl =  paste0(rep(refMotif_true,2),collapse=""), 
                       refMotif_triple = paste0(rep(refMotif_true,3),collapse=""))
correctRef$refMotif_true = as.character(correctRef$refMotif_true)
correctRef_merge = data.frame(refMotif_true = c(correctRef$refMotif_true,correctRef$refMotif_true) ,
                              refMotif_error = c(correctRef$refMotif_dbl,correctRef$refMotif_triple) )
all_fixed = merge(all, correctRef_merge,by.x="refMotif",by.y="refMotif_error",all.x=T)

all_fixed$refMotif.fixed = as.character(all_fixed$refMotif_true);
all_fixed$refMotif.fixed[is.na(all_fixed$refMotif_true)] = as.character(all_fixed$refMotif[is.na(all_fixed$refMotif_true)])

all_fixed$refMotifReverseComplement = stringi::stri_reverse(chartr("ATGC","TACG", all_fixed$refMotif.fixed))

# remove homopolymers and add two "perfect" motif strings where the motif (not 
# reference motif) is replicated to a length of 145 or READ_LEN basepairs total
all_save = all_fixed %>% group_by(refMotifReverseComplement) %>% 
  mutate(lenMotif = length(unique(strsplit(as.character(refMotif.fixed[1]), "")[[1]]))) %>% 
  filter(lenMotif > 1) %>% rowwise() %>% 
  mutate(perfectrep=substr(paste0(rep(Motif,145),collapse=""),1,145),
         perfectrep_long =substr(paste0(rep(Motif,READ_LEN),collapse=""),1,READ_LEN) )

all_save$relevantMotif = ifelse(all_save$revComp == 1,
                                as.character(all_save$refMotifReverseComplement), 
                                as.character(all_save$refMotif.fixed))

all_save_write = all_save[,c("Motif","refMotif.fixed","refMotifReverseComplement","relevantMotif","perfectrep","perfectrep_long")]
colnames(all_save_write) = c("Motif","refMotif","refMotifReverseComplement","relevantMotif","perfectrep","perfectrep_long")
write.table(all_save_write,"motif_2to6bp_correspondence_perfectIRR_noRepeat.txt",row.names=F,col.names=T,quote=F)
