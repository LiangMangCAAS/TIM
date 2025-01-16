# TIM
Transcriptomic Imputation Machine (TIM), an innovative stacking ensemble learning framework designed to improve the accuracy of gene expression imputation and subsequently improve the statistical power of transcriptome-wide association studies (TWAS).


For an example, please refer  to timRun.py.
The example data includes the liver tissue dataset from Chicken-GTEx. The files contain:

Genotype binary files:
    genotype.bed
    genotype.bim
    genotype.fam
    
Annotation file:
    gene_annotation.txt
    
Gene expression data (TPM):
    adjusted_tpm.txt


Note: The example files do not include phenotype data files. For actual use, please prepare the files in the following format.
phenotype.txt

id trait_name
id00001 2.4
id00002 3.7
.
.
.
id001000 2.5








