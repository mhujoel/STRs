set -euo pipefail

CRAM_FILE=$1; READ_FILE=$2
ID=$( basename "$CRAM_FILE" | sed 's/.final.cram//g' );
gunzip ${READ_FILE}
gunzip polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt.gz;

# drop duplicate reads; require SAMFLAG < 256
awk '$4 < 256{print $3}' ${ID}.reads.txt | sort | uniq -c > QNAMES_${ID}.txt

# restrict to potential (non duplicate) IRRs for which the mate is (1) not a potential IRR
# and (2) maps to a locus that is in the polymorphic STR table 
awk 'ARGIND==1{keep[$1] = 1} \
	ARGIND==2 && $1==1{readANC[$2]=1} \
	ARGIND==3 && (FNR == 1 || (keep[$8":"int($9/1e+6 * 10 + 0.5)/10] == 1) && readANC[$3]==1)' \
	polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs_loci.txt \
	QNAMES_${ID}.txt \
	${ID}.reads.txt > reads_sub_${ID}.txt

# restrict to pairs of potential (non duplicate) IRRs; output read name and sequence
awk 'ARGIND==1 && $1==2{readANC[$2]=1} ARGIND==2 && (FNR==1 || readANC[$3]==1){print $3,$10}' \
	QNAMES_${ID}.txt ${ID}.reads.txt > reads_pair_sub_${ID}.txt

# restrict the potential STR table to loci of interest for this individual (locus is 
# potentially represented)
awk 'ARGIND==1{ keep[$8":"int($9/1e+6 * 10 + 0.5)/10] = 1} \
	ARGIND==2 && (FNR == 1 || keep[$3":"int($4/1e+6 * 10 + 0.5)/10] == 1)' \
	reads_sub_${ID}.txt \
	polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt \
	> polyEUR_2to6bp_noHomopolymer_STRs_${ID}.txt

Rscript initial_clean.R ${ID}

# regions to extract to assess if mate is well-anchored:
REGIONS=$(awk 'NR>1{print $3}' initial_clean_${ID}.txt | sort -u | tr '\n' ' ')

# extract regions:
samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
	 ${CRAM_FILE} -b -o ${ID}_mates.bam ${REGIONS}

# restrict to reads that we need to assess quality of:
awk 'ARGIND==1{keep[$5":"$6"_"$4]=1} \
	ARGIND==2 && keep[$3":"$4"_"$1]{print $1,$3,$4,$5}' \
	initial_clean_${ID}.txt \
	<(samtools view ${ID}_mates.bam) > ${ID}_mates_quality.txt

Rscript final_clean.R ${ID}

rm QNAMES_${ID}.txt; rm reads_sub_${ID}.txt;  rm reads_pair_sub_${ID}.txt; 
rm polyEUR_2to6bp_noHomopolymer_STRs_${ID}.txt; rm initial_clean_${ID}.txt; 
rm ${ID}_mates.bam; rm ${ID}_mates_quality.txt;
rm *.crai

gzip ${ID}.reads.txt; 
gzip polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt;
