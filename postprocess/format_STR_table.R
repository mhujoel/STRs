#!/usr/bin/env Rscript
options(echo=FALSE);
param <- commandArgs(trailingOnly=T)

suppressWarnings(suppressMessages(library(dplyr))); options(dplyr.summarise.inform = FALSE)

table1= read.table("allREPEAT_info_simplified_hetEUR.txt",h=T); 
rotations = data.table::fread("motif_2to6bp_correspondence_perfectIRR_noRepeat.txt.gz",h=T)
motif = data.table::fread("motif.txt",h=T)

# add heterozygosity to STR table 
poly = merge(motif,table1[,c("ID","het_EUR")],by="ID")
poly = merge(poly, 
             rotations[,c("Motif", "refMotif", "refMotifReverseComplement")],
             by="Motif")
# add locus to table:
poly$locus_int = paste0(poly$Chrom,":",round(poly$Start/1e+6,digits=1),"_",poly$refMotif); 
poly$repCt = stringr::str_length(poly$refMotif)

# manually add disease relevant STRs that are missing:
df = data.frame(Gene = c("CNBP","COMP","FGF14","NOTCH2NLC","NUTM2BAS1","SAMD12",
                         "STARD7","TBP","XYLT1","BEAN1","DAB1","FXN","MARCHF6",
                         "RAPGEF2","RFC1AAGGG","RFC1ACAGG","TNRC6A","YEATS2"),
                ID=c("chr3:129172576-129172656","chr19:18786034-18786049",
                     "chr13:102161576-102161726","chr1:149390802-149390829",
                     "chr10:79826377-79826403","chr8:118366815-118366913",
                     "chr2:96197066-96197121","chr6:170561907-170562015",
                     "chr16:17470907-17470922","chr16:66490398-66490466",
                     "chr1:57367043-57367118","chr9:69037286-69037304",
                     "chr5:10356347-10356407","chr4:159342526-159342616",
                     "chr4:39348424-39348485","chr4:39348424-39348485",
                     "chr16:24613439-24613529","chr3:183712187-183712222"),
                Motif=c("CAGG","GTC","GAA","GGC","CGG","TGAAA","AAATG","CAG",
                        "GCC","TGGAA","GAAAT","GAA","TTTCA","TTTCA","AAGGG",
                        "ACAGG","TTTCA","TTTCA"),
                Source = "STRipy", het_EUR=NA)

df = merge(df, 
           rotations[,c("Motif", "refMotif", "refMotifReverseComplement")],by="Motif")
df$Motif = as.character(df$Motif)
df$Chrom = stringr::str_split_fixed(df$ID,":",2)[,1];
df$Start = as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(df$ID,":",2)[,2],"-",2)[,1])
df$End = as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(df$ID,":",2)[,2],"-",2)[,2])
df$locus_int = paste0(df$Chrom,":",round(df$Start/1e+6,digits=1),"_", df$refMotif); 
df$repCt = stringr::str_length(df$refMotif)
# -----  add to existing set of polymorphic STR:
ALL = rbind(poly, 
            df[,c("locus_int","ID", "Chrom", "Start", "End","Motif", "refMotif", "refMotifReverseComplement","repCt", "Gene", "Source","het_EUR")])
# ----- write table of STRs:
write.table(ALL,"polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt",row.names=F,col.names=T,quote=F,sep="\t")
