set -euo pipefail
READ_LENGTH=$1;

# Population-specific per locus statistics on allele frequency, heterozygosity, and the
# number of called samples can be found here. Statistics are computed using statSTR from 
# the TRTools package.
wget https://ensemble-tr.s3.us-east-2.amazonaws.com/tables/afreq_tables.zip

# Per locus summary statistics can be downloaded from here. Each table has information on 
# coordinates, repeat unit sequence, and potential overlap with genes listed in GENCODE 
# v22 for repeats in EnsembleTR catalog.
wget https://ensemble-tr.s3.us-east-2.amazonaws.com/tables/repeat_tables.zip

unzip afreq_tables.zip; unzip repeat_tables.zip;

# repeat_info_chr*.csv has columns:
# Chrom Start End Motif ID Gene Source

# afreq_het_chr*.csv has columns:
# Chrom Start ID afreq_AFR het_AFR numcalled_AFR afreq_AMR het_AMR numcalled_AMR afreq_EAS
# het_EAS numcalled_EAS afreq_SAS het_SAS numcalled_SAS afreq_EUR het_EUR numcalled_EUR
#  afreq_H3Africa het_H3Africa numcalled_H3Africa

# restrict to STRs with a repeat motif >= 2 and <= 6 in length and het_EUR > 0
CHR=1
awk -F '\t' 'ARGIND==1 && length($4) >= 2  && length($4) <= 6 {keep[$5] = 1} 
	ARGIND==2 && $17 > 0 &&  (FNR==1 || keep[$3] == 1){print $1,$2,$3,$17}' \
	repeat_info_chr${CHR}.csv afreq_het_chr${CHR}.csv > allREPEAT_info_simplified_hetEUR.txt

for CHR in {2..22}
do
echo $CHR
awk -F '\t' 'ARGIND==1 && length($4) >= 2 && length($4) <= 6 {keep[$5] = 1} 
	ARGIND==2 && $17 > 0 && keep[$3] == 1{print $1,$2,$3,$17}' \
	repeat_info_chr${CHR}.csv afreq_het_chr${CHR}.csv >> allREPEAT_info_simplified_hetEUR.txt
done

# save all 2-6bp motif STR (to add motif)
awk 'NR==1 || (length($4) >= 2  && length($4) <= 6 && $4 != "Motif")' \
	<(cat repeat_info_chr*.csv) > motif.txt

# clean up workspace:
rm repeat_tables.zip; rm afreq_tables.zip; rm repeat_info_chr*; rm afreq_het_chr*

# ----- now in R, first generate a table of containing reference motifs (and reverse 
# complements and rotations thereof) as well as "perfect" IRR sequence.
# requires library: dplyr (assumes it is installed) and installs gtools
Rscript make_motif_table.R $READ_LENGTH
gzip motif_2to6bp_correspondence_perfectIRR_noRepeat.txt 

# Now format the tables to output polymorphic STRs
# requires library: dplyr (assumes it is installed)
Rscript format_STR_table.R
gzip polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt

# extract all loci for which a STR exists (helps filter when assigning IRRs)
awk 'NR>1{print $3":"int($4/1e+6 * 10 + 0.5)/10}' \
	<(zcat polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs.txt.gz) | sort -u  \
	> polyEUR_2to6bp_noHomopolymer_nonRepetitive_DiseaseSTRadded_STRs_loci.txt
	
# final workspace clean up
rm allREPEAT_info_simplified_hetEUR.txt; rm motif.txt