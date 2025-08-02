set -euo pipefail

RAP_DATASET=$1; UKB_FILE=$2; LOCUS=$3; REP=$4; CHR=$5; BP_HG19=$6; IBS_NBRS=$7; 

THREADS=$(nproc)
chmod +x computeIBSpbwt; chmod +x imputeIRRs;
GENETIC_MAP_FILE=genetic_map_hg19_withX.txt.gz

# -----
awk -v l=$LOCUS -v r=$REP '$3==l && $7==r' <(zcat ${UKB_FILE}batch{0..204}_output.txt.gz) > temp.txt
 
# extract relevant covariates: age at recruitment, sex, 10 PCs
# participant.p22001      Genetic sex
# participant.p21022      Age at recruitment
# participant.p22006      Genetic ethnic grouping
# participant.p22009_a{1..10}   Genetic principal components | Array 1..10

dx extract_dataset $RAP_DATASET \
	--fields participant.eid,participant.p22001,participant.p21022,participant.p22006,participant.p22009_a1,participant.p22009_a2,participant.p22009_a3,participant.p22009_a4,participant.p22009_a5,participant.p22009_a6,participant.p22009_a7,participant.p22009_a8,participant.p22009_a9,participant.p22009_a10 \
	-o covariates.csv
		
# could also normalize by coverage (see extract coverage.sh for a way to do this in UKB)
		
# extract relevant SNPs: chr2:47416318:G:A (MSH2) and chr5:80638411:G:T (MSH3)
SNP_DIR="/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)/"
bgenix -g "${SNP_DIR}"ukb21007_c2_b0_v1.bgen  -incl-range chr2:47416318-47416318 > chr2.bgen;
bgenix -g "${SNP_DIR}"ukb21007_c5_b0_v1.bgen  -incl-range chr5:80638411-80638411 > chr5.bgen;

plink2 --bgen chr2.bgen ref-first --sample "${SNP_DIR}"ukb21007_c2_b0_v1.sample \
	--export A --set-missing-var-ids @:#\$r:\$a --new-id-max-allele-len 1000 missing \
	--out chr2

plink2 --bgen chr5.bgen ref-first --sample "${SNP_DIR}"ukb21007_c5_b0_v1.sample \
	--export A --set-missing-var-ids @:#\$r:\$a --new-id-max-allele-len 1000 missing \
	--out chr5

rm chr{2,5}.{bgen,log}

# Compute haplotype neighbors (access to UK Biobank data on the Research Analysis Platform needed) 
# This tool (computeIBSpbwt) outputs the same number of haplotype matches per 
# individual (specified by an input parameter); 

FLANK_BP=20000000;  OUT_FILE=hap_neighbors.temp.txt.gz

# bgen and sample file provided by UKB:
BGEN_FILE=/mnt/project/Bulk/Imputation/Haplotypes/ukb22438_c${CHR}_b0_v2.bgen
SAMPLE_FILE="/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${CHR}_b0_v3.sample"

BGEN_FLANKS=$CHR.$BP_HG19.flanks.bgen
BP_START=$[BP_HG19-FLANK_BP]
if [ $BP_START -lt 0 ]; then BP_START=1; fi

bgenix -g $BGEN_FILE -incl-range :$BP_START-$[BP_HG19+FLANK_BP] > $BGEN_FLANKS

./computeIBSpbwt \
    $CHR $BP_HG19 \
    $BGEN_FLANKS "$SAMPLE_FILE" \
    $GENETIC_MAP_FILE \
    $IBS_NBRS \
    $THREADS \
    $OUT_FILE

rm $BGEN_FLANKS
# -------- 
Rscript make_irr.R

./imputeIRRs temp_h3.txt ${OUT_FILE} IRRs_phased_imp.${LOCUS}.${REP}.h3.txt.gz
./imputeIRRs temp_h4.txt ${OUT_FILE} IRRs_phased_imp.${LOCUS}.${REP}.h4.txt.gz
# ------- 
Rscript compute_associations.R ${LOCUS} ${REP}
gzip phenotypes.*.UKB.txt

rm temp.txt; rm IRRs_phased_imp.${LOCUS}.${REP}.h{3,4}.txt.gz
rm ${OUT_FILE}; rm temp_h3.txt; rm temp_h4.txt; 
rm chr{2,5}.raw; rm covariates.csv;



