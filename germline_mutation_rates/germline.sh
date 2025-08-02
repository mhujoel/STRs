set -euo pipefail

GENOTYPE_FILE=$1; LOCUS=$2; CHR=$3; BP_HG19=$4; IBS_NBRS=$5; ADJ_FACTOR=$6
THREADS=$(nproc)
chmod +x computeIBSpbwt

# restrict to confident (at least 5 spanning reads) heterozygotes with a spread in allele
# lengths >= 3
awk 'NR==1 || (sqrt(($2 - $4)^2) >= 3 && $3 >= 5 && $5 >= 5)' \
	$GENOTYPE_FILE > ${LOCUS}_conf_het_spread.txt
	
# Compute haplotype neighbors (access to UK Biobank data on the Research Analysis Platform needed) 
# This tool (computeIBSpbwt) outputs the same number of haplotype matches per 
# individual (specified by an input parameter); 

FLANK_BP=20000000; 
OUT_FILE=hap_neighbors.${LOCUS}.txt.gz

# bgen and sample file provided by UKB:
BGEN_FILE=/mnt/project/Bulk/Imputation/Haplotypes/ukb22438_c${CHR}_b0_v2.bgen
SAMPLE_FILE="/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${CHR}_b0_v3.sample"

GENETIC_MAP_FILE=genetic_map_hg19_withX.txt.gz

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
# ------- phase STR using haplotype neighbors:
Rscript phase_STR.R ${LOCUS}

# ---- determine top haplotype pairs + ancestral allele
Rscript IBDpair_AncestralAllele.R ${LOCUS}

# ---- compute rates
Rscript compute_GermlineMut_Rates.R ${LOCUS} ${ADJ_FACTOR}

# --- clean up:
rm ${LOCUS}_conf_het_spread.txt phased_${LOCUS}.txt hap_neighbors.${LOCUS}.txt.gz

