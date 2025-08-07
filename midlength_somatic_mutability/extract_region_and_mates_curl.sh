set -euo pipefail

ID=$1
CRAM_URL=$2
CHR=$3
BP_START=$4
BP_END=$5
OUT_DIR=$6
REF_FASTA=$7
GAWK_EXE=awk

mkdir -p $OUT_DIR/cram_slices/
mkdir -p $OUT_DIR/reads/

RETRIES=3; MAX_FILESIZE=40000000 # curl parameters
BP_FLANK=2000 # look for anchors within 2kb of target region
BP_FAR=10000 # threshold for mate to be considered far away if on same chromosome
READ_LEN=151

CRAM_GSURL=$1; 

# AoU:
# ID=$(basename $CRAM_GSURL)
# ID=${ID/.cram/}
# ID=${ID/wgs_/}
# CRAM_URL=$(echo $CRAM_GSURL | sed 's&gs://fc-aou-datasets-controlled&https://fc-aou-datasets-controlled.storage-download.googleapis.com&')
CRAI_URL=$CRAM_URL.crai

# wrapper function that supplies curl arguments to use throughout
function curl_wrapper() {
    curl --max-filesize $MAX_FILESIZE --retry $RETRIES -X GET \
	$@ 2> /dev/null    
        # AoU:
	# -H "Authorization: Bearer $GCS_OAUTH_TOKEN" \
	# -H "X-Goog-User-Project: $GCS_REQUESTER_PAYS_PROJECT" \
}

function echo_cram_footer() {
    echo 0f 00 00 00 ff ff ff ff 0f e0 45 4f 46 00 00 00 00 01 00 05 bd d9 4f 00 01 00 06 06 01 00 01 00 01 00 ee 63 01 4b | xxd -r -p
}

#Download the cram index
curl_wrapper "$CRAI_URL" > $OUT_DIR/cram_slices/$ID.crai

#Get the byte range for the cram header
HEADER_END=$( $GAWK_EXE '{print $4-1; exit}' <(zcat $OUT_DIR/cram_slices/$ID.crai) )

#Download the header
curl_wrapper -r 0-$HEADER_END "$CRAM_URL" > $OUT_DIR/cram_slices/$ID.cram_header

#Look up the ID of the relevant chromosome
TID=$(samtools view -H $OUT_DIR/cram_slices/$ID.cram_header | $GAWK_EXE -v chrom=$CHR '$1=="@SQ" {tid1++; if ($2=="SN:"chrom) print tid1-1;}')

#Determine the start and end positions of blocks overlapping the region of interest
# NOTE: this code doesn't handle the corner case in which the cram slice spans multiple contigs
#       see cram_slice.sh for a robust version
read -r RANGE_START RANGE_END < \
    <(zcat $OUT_DIR/cram_slices/$ID.crai \
    | $GAWK_EXE -v tid=$TID -v bpstart=$[BP_START-BP_FLANK] -v bpend=$[BP_END+BP_FLANK] '
{ # process cram slice
  if (start_printed==1 && (($1==tid && $2>bpend) || $1!=tid)) { # at first slice after region
    print $4-1;
    exit;
  }
  if (start_printed!=1 && $1==tid && $2+$3>bpstart) { # found first slice that intersects region
    printf $4" ";
    start_printed=1;
  }
}')

#Download the blocks
curl_wrapper -r $RANGE_START-$RANGE_END $CRAM_URL > $OUT_DIR/cram_slices/$ID.cram_blocks

#Concatenate the header, blocks, and footer
CRAM_SLICE=$OUT_DIR/cram_slices/$ID.slice.cram
cat $OUT_DIR/cram_slices/$ID.{cram_header,cram_blocks} > $CRAM_SLICE
echo_cram_footer >> $CRAM_SLICE

#Find reads with no bad flags (but keeping supplementary alignments) that satisfy:
#- overlap repeat
#- mates aligned far away
touch $OUT_DIR/reads/$ID.overlapping_reads.txt
touch $OUT_DIR/reads/$ID.anchors_with_far_mates.txt
samtools view -T $REF_FASTA -F 0x700 -F 0x4 -F 0x8 $CRAM_SLICE \
    | awk -v bpstart=$BP_START -v bpend=$BP_END -v bpflank=$BP_FLANK -v bpfar=$BP_FAR -v readlen=$READ_LEN -v id=$ID -v out_dir=$OUT_DIR '
$4 >= bpstart-readlen+1 && $4 <= bpend && $2 < 2048 {
  print | "cut -f-11 > "out_dir"/reads/"id".overlapping_reads.txt"
}
$7 != "=" || $8 < bpstart-bpfar-readlen+1 || $8 > bpend+bpfar {
  if (($4 > bpstart-bpflank-readlen+1 && $4 < bpend+bpflank && $2 < 2048 && $5 >= 30) || ($4 > bpstart-readlen+1 && $4 < bpend))
    print | "cut -f-11 > "out_dir"/reads/"id".anchors_with_far_mates.txt"
}'

function extract_mates() {
    READS_FILE=$1
    OUTPUT_FILE=$2

    #Use awk to identify additional cram slices that contain mates of reads; build mate cram file
    cp $OUT_DIR/cram_slices/$ID.cram_header $OUT_DIR/cram_slices/$ID.mates.cram
    $GAWK_EXE -v chr=$CHR -v readlen=$READ_LEN '
ARGIND==1 && $1=="@SQ" { # process cram header: create lookup table (contig name -> tid)
  tid1++;
  contig2tid[substr($2,4)] = tid1-1;
}
ARGIND==2 { # process anchors with far mates: record positions of mates
  mate_chr = $2; if (mate_chr=="=") mate_chr = $1;
  mate_tid = contig2tid[mate_chr]; if (mate_tid=="") mate_tid = -1;
  mate_pos = $3;
  bpTargetCtPerTid[mate_tid]++;
  bpTargetsByTid[mate_tid,bpTargetCtPerTid[mate_tid]] = mate_pos;
}
ARGIND==3 { # process cram index; write byte ranges for records overlapping the above
  if (first_byte) { # previous line contained a pos of interest => print last byte of slice
    if ($4 == first_byte) # ... unless this record is still in the same slice; skip to next if so
      next;
    print $4-1;
  }
  tid = $1;
  slice_start = $2;
  slice_end = $2+$3-1;
  first_byte = 0; # set to nonzero below if this slice is to be extracted
  for (i = 1; i <= bpTargetCtPerTid[tid]; i++)
    if (slice_start <= bpTargetsByTid[tid,i] && bpTargetsByTid[tid,i] <= slice_end)
      first_byte = $4;
  if (slice_end - slice_start < 2*readlen)
    first_byte = 0; # skip slice if coverage is excessively high
  if (first_byte)
    printf first_byte" "; # print first byte of this slice
}' \
    <(samtools view -H $OUT_DIR/cram_slices/$ID.cram_header) \
    <(cut -f3,7-8 $READS_FILE | sort | uniq) \
    <(zcat $OUT_DIR/cram_slices/$ID.crai) \
    | while read -r RANGE_START RANGE_END
    do
	curl_wrapper -r $RANGE_START-$RANGE_END $CRAM_URL >> $OUT_DIR/cram_slices/$ID.mates.cram
    done
    echo_cram_footer >> $OUT_DIR/cram_slices/$ID.mates.cram

    #Extract mates with read names matching query set
    samtools view -T $REF_FASTA -F 0xF00 -N <(cut -f1 $READS_FILE) \
	$OUT_DIR/cram_slices/$ID.mates.cram \
	| cut -f-11 > $OUTPUT_FILE
}

#Extract mates of anchors with far mates
extract_mates $OUT_DIR/reads/$ID.anchors_with_far_mates.txt $OUT_DIR/reads/$ID.far_mates_or_primary.txt

#Check whether any mates of far mates haven't yet been extracted
#(this happens if the original read was a supplementary alignment, and the primary alignment wasn't picked up in the search for far mates)
awk 'ARGIND==1 && $2 >= 2048 { read_was_supp[$1]=1 }
ARGIND==2 { read_obs_count[$1]++; read_info[$1]=$0; }
END {
  for (read in read_info)
    if (read_obs_count[$1]==1 && read_was_supp[$1])
      print read_info[$1];
}' \
    $OUT_DIR/reads/$ID.anchors_with_far_mates.txt \
    $OUT_DIR/reads/$ID.far_mates_or_primary.txt \
    > $OUT_DIR/reads/$ID.far_mates_missing_primary.txt

#Extract still-at-large mates of far mates of supplementary alignments (if needed)
if [ -s $OUT_DIR/reads/$ID.far_mates_missing_primary.txt ]; then
    extract_mates $OUT_DIR/reads/$ID.far_mates_missing_primary.txt $OUT_DIR/reads/$ID.supp_with_far_mates_primary.txt
    $GAWK_EXE -v id=$ID '
ARGIND==1 {nFarMatesMissingPrimary++}
ARGIND==2 {nPrimaryRetrieved++}
END { printf("%s: retrieved %d off %d missing primary alignments not found in far mates\n",id,nFarMatesMissingPrimary,nPrimaryRetrieved); }' \
        $OUT_DIR/reads/$ID.far_mates_missing_primary.txt \
        $OUT_DIR/reads/$ID.supp_with_far_mates_primary.txt    

    #Append retrieved primary alignments to overlapping reads file
    cat $OUT_DIR/reads/$ID.supp_with_far_mates_primary.txt >> $OUT_DIR/reads/$ID.overlapping_reads.txt
    rm $OUT_DIR/reads/$ID.supp_with_far_mates_primary.txt
fi

#Summarize results
$GAWK_EXE -v id=$ID '
ARGIND==1 {nOverlap++}
ARGIND==2 {nAnchors++; anchorChrPos[$1]=$7" "$8}
ARGIND==3 {nFarMates++; farMateFound[$1]=1}
END {
  for (read in anchorChrPos)
    if (!farMateFound[read])
      printf("%s: skipped mate at %s\n",id,anchorChrPos[read]);
  printf("%s: extracted %d overlapping reads, %d of %d far mates\n",id,nOverlap,nFarMates,nAnchors);
}' \
    $OUT_DIR/reads/$ID.overlapping_reads.txt \
    $OUT_DIR/reads/$ID.anchors_with_far_mates.txt \
    <(cut -f1 $OUT_DIR/reads/$ID.far_mates_or_primary.txt | sort | uniq) \
    | sort

#Clean up
rm $OUT_DIR/cram_slices/$ID.{cram_blocks,cram_header,crai,slice.cram,mates.cram} \
    $OUT_DIR/reads/$ID.far_mates_missing_primary.txt
