To extract information from confident heterozygotes necessary to compute rates of somatic expansion of short alleles two scripts are needed ``somatic.sh`` and ``short_somatic_perIndividual.R``. 

Then in a directory with a reference genome fasta file and the 2 necessary scripts, one can run the script with the following command:
```
bash somatic.sh PATH_TO_CRAM_FILE LOCUS_CONSIDERED
```
where ``PATH_TO_CRAM_FILE`` is the path to a cram file of interest, and ``LOCUS_CONSIDERED`` is a string of the following format: ``LOCUS_CHR_STARTBP_ENDBP_STARTstr_ENDstr_REP_L6,9_R,`` where ``LOCUS, CHR, STARTBP, ENDBP`` are the locus name and location; ``STARTstr`` and ``ENDstr`` are sequences on the flank of the repeat (6 base pairs + 3 base pair of repeat sequence); ``REP`` is the repeat motif; ``LFLANK`` and ``RFLANK`` are of the format: ``{L,R}[comma-separated list of "key" flank bases]`` or ``{L,R},`` if no "key" flank bases; where a key base is one where a base error could result in there appearing to be one more or one fewer repeat unit  

For example, for TCF4 one could run:
```
bash somatic.sh ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240216/HG00259.final.cram TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,
```
The output is ``summary_somatic_HG00259.txt`` (included). 

The columns are:

``V1`` : ID of individual

``V3`` : length of repeat sequence (in base pairs)

``V4`` : repeat sequence of read

``jump`` : jump (in repeat units) of sequence from main allele (-2,-1,1,2); for main alleles jump is coded as NA

``fromLEN`` : length of repeat sequence (in base pairs) from which the read may have potentially somatically expanded/contracted from

``A1len A1midconsensus`` : repeat length (in base pairs) and consensus sequence of main allele 1

``A2len A2midconsensus`` : repeat length (in base pairs) and consensus sequence of main allele 2

``somatic_neg2 somatic_neg1 somatic_pos1 somatic_pos2``:  indication if the read passes the "potentially-somatic" filters (1 if passes,0 otherwise) if the read originated from a -2, -1, +1, +2 jump 


Note: scripts are currently written for 1000 Genomes formatted cram paths; for example ``ID=$( basename "$FILE" | sed 's/.final.cram//g' )``; for UK Biobank cram files this would need to be modified. 

How to compile and use these ``summary_somatic_*.txt`` and what you would get if you ran this across all of UKB is included below. We assume the compiled ``summary_somatic_*.txt`` files for UKB are concatenated within ``all.txt``
For each allele length + jump -- we need 
(1) number of somatic reads 
(2) number of individuals considered [n per allele] 

```
tab_sum  = read.table("all.txt",h=T)

library(dplyr)
sem <- function(x) {
  return(sd(x)/sqrt(length(x)))
}

repLen = 3
# number of individuals per allele:
number_per_allele = unique(tab_sum[,c("V1","A1len","A1midconsensus","A2len","A2midconsensus")])
allele_cts = data.frame(A = c(number_per_allele$A1len,number_per_allele$A2len),Aseq = c(as.character(number_per_allele$A1midconsensus),as.character(number_per_allele$A2midconsensus))) %>% group_by(A, Aseq) %>% summarize(n=n())

# somatic read counts: get the number of reads that appear likely a somatic expansion (+1,+2) or contraction (-1,-2) from each allele length / repeat sequence

somatic_cts = tab_sum %>% filter(!is.na(jump) & !is.na(somatic_neg2) & !is.na(somatic_neg1) & !is.na(somatic_pos1) & !is.na(somatic_pos2)) %>%  rowwise() %>% 
  mutate(relevantSOM = ifelse(jump==-2,somatic_neg2,ifelse(jump==-1,somatic_neg1,ifelse(jump==1,somatic_pos1,somatic_pos2))),
         relevantALLELE = ifelse(fromLEN==A1len,as.character(A1midconsensus),as.character(A2midconsensus))) %>%
  group_by(jump,fromLEN,relevantALLELE) %>%  summarize(nSomatic_reads = sum(relevantSOM))

# number of reads expected : number of trustworthy reads (from reads originating from main alleles) that passed likely-somatic filters 
# e.g. for +1 expansion from an 18 allele (to a 19 allele) we compute the mean number of main allele 19 reads that pass likely-somatic filters as if it was a +1 jump

denominator_cts = tab_sum %>% filter(is.na(jump) & !is.na(somatic_neg2) & !is.na(somatic_neg1) & !is.na(somatic_pos1) & !is.na(somatic_pos2)) %>%   group_by(V1,V3) %>% 
  summarize(jump=c(-2,-1,1,2), nREADs = c(sum(somatic_neg2),sum(somatic_neg1),sum(somatic_pos1),sum(somatic_pos2))) %>%
  group_by(V3,jump) %>% summarize(nDENOMINATOR=n(),DENOMINATOR = mean(nREADs),SEM_DENOMINATOR = sem(nREADs)) %>%
  rowwise() %>% mutate(fromLEN = V3 - jump*repLen)

df_complete = merge(somatic_cts, denominator_cts[,c("jump","fromLEN","nDENOMINATOR", "DENOMINATOR", "SEM_DENOMINATOR")], by=c("jump","fromLEN"))
df_complete = merge(df_complete, allele_cts, by.x=c("fromLEN","relevantALLELE"),by.y=c("A","Aseq"))

rates = df_complete %>% group_by(jump,fromLEN,relevantALLELE,nDENOMINATOR,n) %>% summarize(rate = (nSomatic_reads/n)/DENOMINATOR)

as.data.frame(rates[rates$jump == 1,] %>% arrange(fromLEN) %>% filter(n > 500) %>% rowwise() %>% mutate(Alen=fromLEN/3,relevantA = stringr::str_replace_all(relevantALLELE,"AGC",".")) %>% ungroup() %>% select(c(fromLEN,Alen,relevantA,rate)))

   fromLEN Alen                            relevantA         rate
1       36   12                         ............ 0.0005799788
2       39   13                        ............. 0.0004785547
3       42   14                       .............. 0.0014116528
4       45   15                      ............... 0.0028429565
5       48   16                     ................ 0.0034989797
6       51   17                    ................. 0.0034487537
7       54   18                 ........ACC......... 0.0003915314
8       54   18                   .................. 0.0053746673
9       57   19                  ................... 0.0054112537
10      57   19                ...AGG............... 0.0012762308
11      60   20                 .................... 0.0073812076
12      63   21                ..................... 0.0081248854
13      66   22               ...................... 0.0089533828
14      69   23              ....................... 0.0111291469
15      72   24             ........................ 0.0107681726
16      75   25            ......................... 0.0114192670
17      78   26           .......................... 0.0122385328
18      81   27          ........................... 0.0138454948
19      84   28         ............................ 0.0119092871
20      87   29        ............................. 0.0137258301
21      90   30       .............................. 0.0124575164
22      93   31      ............................... 0.0176779860
23      96   32     ................................ 0.0145977528
24      99   33    ................................. 0.0144791165
25     102   34   .................................. 0.0131540600
26     105   35  ................................... 0.0151510166
27     108   36 .................................... 0.0000000000
```
