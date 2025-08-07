To determine which loci have signal of somatic instability in UK Biobank, this script, given a ``LOCUS`` to test, will:

1) extract IRRs of given locus across UK Biobank; assumes the IRRs have been extracted and assigned as is described in the ``postprocess`` folder and collapsed across individuals with no header line. Currently the script assumes these files have been batched into 205 batches, labeled 0 to 204, and gzipped in the directory: ``${UKB_FILE}``. I.e. the script extracts individuals with IRRs at the locus from the files: ``${UKB_FILE}batch{0..204}_output.txt.gz``

2) Extracts, from ``$RAP_DATASET``, Genetic sex (``participant.p22001``), Age at recruitment (``participant.p21022``), Genetic ethnic grouping (``participant.p22006``), and 10 genetic principal components (``participant.p22009_a{1..10}``)

3) Extracts chr2:47416318:G:A (MSH2) and chr5:80638411:G:T (MSH3) from "/Bulk/Imputation/Imputation from genotype (TOPmed)/" bgen files

4) For each haplotype, extracts a given number (``$IBS_NBRS``) of IBS matches at the given locus and uses these to impute the number of IRRs at the locus for each individual's haplotype.

** Note the script currently does not adjust the number of IRRs for coverage; to adjust for coverage the scripts have to be modified (notes on how to do so are within the scripts). Some notes on extracting coverage from DRAGEN output is provided in ``extract_coverage.sh``.

5) Residualizes and crops the number of IRRs (with a Hamming distance of 3 or 4) and associates this measurement with age, chr2:47416318:G:A (MSH2) and chr5:80638411:G:T (MSH3) conditioning on 10 genetic principal components and genetic sex.

To determine if a locus has signal of somatic instability in UK Biobank:
```
	bash longSignal.sh $RAP_DATASET $UKB_FILE $LOCUS $REP $CHR $BP_HG19 $IBS_NBRS
```
in a directory with ``genetic_map_hg19_withX.txt.gz``, the executables ``computeIBSpbwt`` and ``imputeIRRs``, and the scripts ``compute_associations.R`` and ``make_irr.R``. 

The outputs are:

``associations.$LOCUS.$REP.UKB.txt`` : associations of residualized and cropped IRR counts with age, chr2:47416318:G:A (MSH2) and chr5:80638411:G:T (MSH3); associations restricted to samples who self-identified as 'White British' according to Field 21000 and have very similar genetic ancestry based on a principal components analysis of the genotypes.

``phenotypes.$LOCUS.$REP.UKB.txt.gz`` : residualized and cropped IRR counts and covariates among samples who self-identified as 'White British' according to Field 21000 and have very similar genetic ancestry based on a principal components analysis of the genotypes. This file can be used for downstream GWAS studies on somatic instability of the given repeat locus.

For example, with ``LOCUS=chr19:14774675-14774689; REP=AAAG; IBS_NBRS=10; CHR=19; BP_HG19=14885487`` running the script results in 
``associations.chr19:14774675-14774689.AAAG.UKB.txt`` (below)

```
                   LOCUS  REP DEF Nnonzero MSH2beta MSH2se    MSH2p MSH3beta  MSH3se MSH3p AGEbeta   AGEse    AGEp
 chr19:14774675-14774689 AAAG  H3   191889     0.56  0.021 3.2e-156     0.31  0.0052     0  0.0054 0.00045 1.5e-33
 chr19:14774675-14774689 AAAG  H4   197909     0.60  0.022 4.3e-163     0.33  0.0055     0  0.0057 0.00047 7.0e-35
```

and ``phenotypes.chr19:14774675-14774689.AAAG.UKB.txt.gz``


