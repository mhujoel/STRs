Given the DRAGEN STR files in UKB all STRs can be extracted using the following 2 scripts (and launch script); where this script extracts:
```
REPCN "Number of repeat units spanned by the allele"
ADSP "Number of spanning reads consistent with the allele"
REPCI "Confidence interval for REPCN"
SO "Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat" 
```
It assumes the scripts are in a directory ``UKB_RAP_PATH_TO_SCRIPTS`` and the output goes to ``${DIR}/all_DRAGEN_extracted/``. 

```
for BATCH in {10..60}
do
echo -n ${BATCH}
dx run swiss-army-knife \
-iin="${UKB_RAP_PATH_TO_SCRIPTS}"/extract_all_dragen_batch.sh -iin="${UKB_RAP_PATH_TO_SCRIPTS}"/extract_all_dragen.sh \
-icmd="set -euo pipefail && BATCH=$BATCH /usr/bin/time -v bash extract_all_dragen_batch.sh 2>&1 | tee batch${BATCH}.log" \
--instance-type mem2_ssd1_v2_x2 \
--destination "${DIR}"/all_DRAGEN_extracted/ \
--priority low --brief \
--name DRAGEN_ALL_$BATCH \
-y
done

# ------------------ extract_all_dragen_batch.sh -------------------- #
# set on input: $BATCH
set -euo pipefail
sudo apt-get --yes install parallel

DRAGEN_DIR="/mnt/project/Bulk/DRAGEN WGS/Whole genome STR call files (DRAGEN) [500k release]/$BATCH/"
OUT=batch${BATCH}_DRAGEN.txt
set +e # OK if some jobs fail; status will appear in job log
ls "$DRAGEN_DIR"/*.gz | parallel --retries 3 --joblog /dev/stderr -j 8 bash extract_all_dragen.sh {} >> $OUT
set -e
gzip ${OUT}
# ------------------ extract_all_dragen.sh -------------------- #
set -euo pipefail
FILE="$1"; ID=$( basename "$FILE" ); ID=${ID:0:7};
bcftools query -f '%VARID\t%RU[\t%REPCN\t%ADSP\t%REPCI\t%SO]\n' "$FILE" | awk -v iid=$ID '{print iid,$0}' | sed 's=/= =g'
# ------------------------------------------------------------------#
```


Genotypes for TCF4 can then be extracted as follows:
```
echo -e 'IID A1 A1rd A2 A2rd'  > TCF4_DRAGEN_GENOTYPES.txt
zcat ${DIR}/all_DRAGEN_extracted/batch*_DRAGEN.txt.gz | awk '$2=="TCF4" && $10 == "SPANNING" && $11 == "SPANNING" {print $1,$4,$6, $5,$7}' >> TCF4_DRAGEN_GENOTYPES.txt
```
