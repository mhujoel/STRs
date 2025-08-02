set -euo pipefail
FILE="$1"; ID=$( basename "$FILE" | sed 's/.final.cram//g' ); DAT=$2;

CHR=$(echo $DAT | cut -d'_' -f 2); START=$(echo $DAT | cut -d'_' -f 3); END=$(echo $DAT | cut -d'_' -f 4);
STARTSEQ=$(echo $DAT | cut -d'_' -f 5); ENDSEQ=$(echo $DAT | cut -d'_' -f 6);

samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa -F 0x4 -F 0x100 -F 0x200 -F 0x400 -F 0x800 "$FILE" ${CHR}:${START}-${END} \
| awk -v s=$STARTSEQ -v e=$ENDSEQ -v iid=$ID '
BEGIN {STARTSEQ_TOADD = substr(s,7,3); ENDSEQ_TOADD=substr(e,1,3)}
$5 >= 30 && $7=="=" {
start=split($10,a,s);end=split($10,b,e); 
if (start == 2 && end == 2) {
split(a[2],segMID,e);
print iid,(and($2,16)==0?"forward":"reverse"),length(STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD), STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD,$10,$11
}
}' > IID_${ID}.txt

nREADS=$(cat IID_${ID}.txt | wc -l)

if [ $nREADS -gt 0 ]; then
Rscript short_somatic_perIndividual.R $ID $DAT
fi

rm IID_${ID}.txt; rm *.crai
