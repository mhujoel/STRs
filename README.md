# Code and scripts for methods described in "Causes and consequences of DNA repeat expansions among 900,000 biobank participants"

We provide the following components:

1)	Standalone tool (``extractLongSTRs``) for efficiently extracting all reads with high repeat content (i.e., potential IRRs) from a WGS cram file generated using any standard read-mapper such as bwa or DRAGEN. A GitHub repository for this tool that includes C++ code, binary executables, documentation, and an example can be found at: XXX.
2)	Pipeline of bash and R scripts that post-processes extracted reads to identify IRRs that match known STR loci (from Zaiei Jam et al. 2023). We have included an example using publicly available 1000 Genomes Project data. See ``postprocess``.
3)	Pipeline for further processing IRR measurements into somatic expansion phenotypes and running GWAS and PheWAS on the UKB-RAP platform. This pipeline includes code that imputes inherited allele lengths (which are used to generate somatic expansion phenotypes) using phased SNP-haplotypes available on UKB-RAP. See ``processIRR_somaticexpansion``
4)	Pipeline for assessing somatic mutability of short alleles at a given STR locus (using base quality information to filter potential PCR stutter errors) including an example using 1000 Genomes Project data. See ``short_somatic_mutability``.
5)	Pipeline for estimating germline mutation rates of alleles at a given STR locus (using analysis of IBD segments) using phased SNP-haplotypes available on UKB-RAP. We have provided an example using DRAGEN STR genotypes available on UKB-RAP. See ``germline_mutation_rates``.
6)	Pipeline for quantifying mean somatic expansion of mid-length alleles of AAAG repeats at chr2:232.4Mb and chr19:14.8Mb including an example using 1000 Genomes Project data. See ``midlength_somatic_mutability``.
7)	C++ files are provided that implement estimation of inherited allele lengths (via imputation, as invoked in the above pipelines) and estimation of mean somatic expansion of mid-length AAAG repeat alleles at 2 loci. See ``cpp_files``.

   
The reference fasta can be downloaded here: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.{fa,fa.fai}

The genetic map file can be downloaded here: https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz 
