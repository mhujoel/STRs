Given a text file of genotypes within a large cohort (shown here for UK Biobank), the script ``germline.sh`` computes germline instability rates. In more detail:
1) It restricts to confident (at least 5 spanning reads) heterozygotes with a spread in allele lengths >= 3 
2) Extracts haplotype matches at a given locus (** currently coded for UK Biobank data)
3) Uses haplotype matches to phase STR genotypes
4) Restricts to an individual's top haplotype match and computes their allele length difference (delta) and the ancestral allele of the pair (if possible)
5) Computes germline mutation rates (adjusting for ascertainment of only the longest IBD match for each haplotype)

Note that an input is the text file of genotypes which contains 5 columns of (1) individual label; IID (2 & 3) allele 1 repeat length (A1) and its spanning read support; and (4 & 5) allele 2 repeat length (A2) and its spanning read support with column labels: ``IID,A1,A1rd,A2,A2rd``. Allele repeat lengths should be in units of repeat units.

To compute germline instability rates in UK Biobank:
```
	bash germline.sh PATH_TO_GENOTYPE_FILE LOCUS CHR BPHG19 NUMIBS ADJ_FACTOR
```
in a directory with ``genetic_map_hg19_withX.txt.gz`` and ``autosome10k_ibdne.ne.txt``, the executable ``computeIBSpbwt``, and the scripts ``phase_STR.R`` ``IBDpair_AncestralAllele.R`` and ``compute_GermlineMut_Rates.R``. ``NUMIBS`` is the number of haplotype matches per individual to extract in whole cohort; note that when phasing and looking at haplotype matches we restrict to individuals that are confident heterozygotes - and thus some haplotype matches may be dropped.

``autosome10k_ibdne.ne.txt`` contain estimates of effective population size in UK Biobank kindly provided by R. Cai and S. Browning (Cai, R., Browning, B. L. & Browning, S. R. Identity-by-descent-based estimation of the X chromosome effective population size with application to sex-specific demographic history. G3 (2023)).

The outputs are:
 
``germline_mut_${LOCUS}_extended_info.txt``  : for each individual haplotype: top haplotype neighbor, cM shared (and flank amount), ancestral allele, number considered for ancestral allele and number that agreed on the determined ancestral allele. 

``germline_mut_${LOCUS}_unique_pairs.txt `` : unique ID1-ID2 pairs with cM > 5, cM flank > 0.5, an ancestral allele determine (at least 3 agree and at least 5 were considered); output is the ID1, ID2, cM, ancestral allele, and delta (-9 if neither in pair matches ancestral allele, otherwise is difference to ancestral allele, if any, of mutated allele)

``germline_mut_rates_${LOCUS}.txt`` : mutation rate of +/- 1,2 jumps per allele with at least 500 individuals. These rates are adjusted for the fact that we ascertained the longest IBS match (in UKB we multiply by 1.204548)

So for example, for TCF4 (hg38; chr18:55586154-55586228; hg19: chr18:53253385-53253459) one could run:
```
bash germline.sh TCF4_DRAGEN_GENOTYPES.txt TCF4 18 53253385 20 1.204548
```
See "README_DRAGEN_STR" for notes on how to extract DRAGEN STR calls in UK Biobank.

We have made available the output of this command that is non-identifiable (the germline mutation rates).

