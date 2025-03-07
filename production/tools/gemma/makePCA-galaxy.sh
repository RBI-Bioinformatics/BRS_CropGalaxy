#!/bin/bash

# This tool calls PLINK to run Principal Component Analysis and
# converts the output into a form accepted as a GEMMA covariate file

# Inputs
#bed=$1  # PLINK BED file
#bim=$2  # PLINK BIM file
#fam=$3  # PLINK FAM file
input=$1

n_pc=${2:-10}  # Number of PCs to return,  suggested default 10  (then extract smaller sets if needed).

# Output basename. This tool creates several files.
outbase=pca

# A function to convert to a GEMMA covariate file. Note that the first column will be always 1 (intercept).
pca2gemma(){
  pcafile=${1?"Please provide an input eigenvec file"}
  skip=${2:-0}
  >&2  echo "Will skip $skip lines"
  awk -v skip="$skip" 'NR > skip { $1=""; $2 = ""; print 1,$0 }' $pcafile | tr -s " "
}

## Main computation
#plink --bed $bed --bim $bim --fam $fam \
plink --bfile $input \
  --keep-allele-order \
  --pca $n_pc header tabs \
  --out $outbase \
  --allow-extra-chr --allow-no-sex

rm *.nosex # remove that file - not needed

# The created  PCA eigenvec file can be used for plink GLM, but for GEMMA covariate we need to transform:
pca2gemma $outbase.eigenvec 1  > pca_for_gemma.txt

# Outputs renaming
mv $outbase.eigenvec PC.txt
mv $outbase.eigenval Eigenvalues.txt

# For user convenience, if the requested number of PC is large, extract commonly used subsets
# This is done because often users try to run GWAS with different numbers of PCs computed on the same dataset.
# It's better to compute many PCs once, then extract smaller numbers of top PCs
if [[ $n_pc -gt 5 ]] ; then
 # cut -d" " -f1-2 pca_for_gemma.txt > PC1.txt
 cut -d" " -f1-3 pca_for_gemma.txt > PC1-2.txt
 # cut -d" " -f1-4 pca_for_gemma.txt > PC1-3.txt
 # cut -d" " -f1-5 pca_for_gemma.txt > PC1-4.txt
 cut -d" " -f1-6 pca_for_gemma.txt > PC1-5.txt
fi
######################
## Outputs:
# PC.txt             (tabular)
# Eigenvalues.txt    (1-colum tabular, no header)
# pca_for_gemma.txt  ( Npc+1 column space separated, txt)
# (possible output:)
# PC1-2.txt          (3-colum txt)
# PC1-5.txt          (6-colum txt)
###################
