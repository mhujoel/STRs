A BGEN file containing genotypes corresponding to if someone carriers an expanded (e.g. at least 1 IRR) repeat at a locus can be constructed in R, first writing bed/bim/fam files using ``genio::write_plink()`` function; this can be done per batch then merged using ``plink``;  sample command to merge across files:
```
for BATCH in {0..204}
do
echo ${DIR}/H3_anchor.batch${BATCH}.{bed,bim,fam} >> all_my_files.txt
done

plink --merge-list all_my_files.txt --make-bed --out H3_anchor_bfiles
```
Then ``plink2`` can be used to create a ``bgen`` file:
```
plink2 --bfile H3_anchor_bfiles --export bgen-1.2 bits=8 --out H3_anchor 
```
NOTE: after creating bgen file, the batch-level bed,bim,fam levels can be removed!

A PheWAS can then be conducted by modifying the following code: 
```
BGEN_PREFIX=$PATH_TO_BGEN
DIR_BOLT=$PATH_TO_BOLT
tar -zxf $DIR_BOLT/BOLT-LMM_v2.5.tar.gz
NUM_PCS=20

./BOLT-LMM_v2.5/bolt \
    --bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c21_b0_v2" \
    --remove=$PATH_TO_REMOVE_FILE
    --phenoFile=$PATH_TO_PHENOTYPE_FILE --phenoCol=${PHENO} \
    --covarFile=$PATH_TO_COVARIATE_FILE \
    --covarCol=cov_SEX_GENETIC --qCovarCol=cov_AGE --qCovarCol=cov_AGE_SQ \
    --qCovarCol=ukbPC{1:$NUM_PCS} \
    --numThreads=$(nproc) \
    --statsFile=/dev/null \
    --bgenFile="$BGEN_PREFIX.bgen" --sampleFile="$BGEN_PREFIX.sample" \
    --bgenMinINFO=0.3 \
    --statsFileBgenSnps=H3_anchor.$PHENO.bgen.stats.gz \
    --verboseStats
```
