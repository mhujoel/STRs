# the following batch script can be run for batches 10..60:
# --------- coverage_batch.sh 
# set on input: $BATCH
set -euo pipefail
sudo apt-get --yes install parallel

SUPP_DIR="/mnt/project/Bulk/DRAGEN WGS/Whole genome supplementary files (DRAGEN) [500k release]/$BATCH"

OUTPUT_COV=COVERAGE_${BATCH}.txt
echo -e 'ID\tAVG' > ${OUTPUT_COV}

set +e # OK if some jobs fail; status will appear in job log
ls "$SUPP_DIR"/*.zip | parallel --joblog /dev/stderr -j 4 bash coverage.sh {} >> ${OUTPUT_COV}
set -e

# -------- coverage.sh
set -euo pipefail
FILE="$1"; ID=$( basename "$FILE" ); ID=${ID:0:7}
unzip -qq "${FILE}" -d temp_${ID}
COV=$(grep "COVERAGE SUMMARY,,Average autosomal coverage over genome," temp_${ID}/*.wgs_coverage_metrics.csv | sed 's/COVERAGE SUMMARY,,Average autosomal coverage over genome,//g')
rm -r temp_${ID}/
echo -e $ID'\t'$COV

	