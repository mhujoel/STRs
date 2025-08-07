set -euo pipefail

REF_FASTA=GRCh38_full_analysis_set_plus_decoy_hla.fa # fa+fai unpacked or symlinked to working dir
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

READ_LEN=150 # for 1000G 30x

for LOCUS_INFO in "chr2 232389717 232389783 chr2_232_4Mb" "chr19 14774675 14774689 chr19_14_8Mb"
do
    read -r CHR BP_START BP_END OUT_PREFIX <<< "$LOCUS_INFO"
    SUMMARY=$OUT_PREFIX.1000G_GBR.allele.mean_somatic_expansion.txt
    rm -f $SUMMARY
    while read -r ID CRAM_URL
    do
	# echo $ID $CRAM_URL
	bash extract_region_and_mates_curl.sh $ID $CRAM_URL $CHR $BP_START $BP_END $OUT_PREFIX $REF_FASTA
	cat $OUT_PREFIX/reads/$ID.{overlapping_reads,far_mates_or_primary}.txt \
	    | ./find_repeats_$OUT_PREFIX $ID $READ_LEN \
	    >> $SUMMARY
    done < 1000G_GBR_cram_URLs.txt
done
