# BGGN239
#### Data and instructions for BGGN 239 in-class exercises

```
git clone https://github.com/TiffanyAmariuta/BGGN239
```

Genes of interest for Day 1 exercises: 
ENSG00000158864.12 for Brain_Cerebellum
ENSG00000116704.7 for Whole_Blood

Diseaes of interest for Day 1 exercises: 
PASS_Alzheimers_deRojas2021.sumstats (Alzheimer's)
PASS_UC_deLange2017.sumstats (Ulcerative Colitis)

GTEx gene expression files: 
dimensions: genes on the rows, people on the columns (deidentified IDs)
genes on the rows are a subset of all genes that would be included in the full GTEx dataset. 
Each file is specific to a different tissue because different people were used for brain samples vs blood samples: 
GTEx_matrix_blood.txt
GTEx_matrix_brain.txt

GTEx covariate files: 
dimensions: people on the rows (matches the order of the GTEx columns) and covariates on the columns (ancestry PCs, PEER factors, sex, and technology of sequencing)
Each file is specific to a different tissue because different people were used for brain samples vs blood samples: 
covar_blood.txt
covar_brain.txt

GTEx plink files: 
GTEx genetic data files (plink) come in a trio (bed/bim/fam) for the entire GTEx dataset. These files contain genotypes for all SNPs for all people in the entire dataset (e.g. across tissues). 
I have subset these plink files into just the SNPs and just the people you'll need for today's exercises. 
Each file is specific to a gene's unique cis-window (subset of nearby SNPs) and the people as found in the covar_* and GTEx_matrix_* files above: 
ENSG00000158864.12_brain (with .bed, .bim, .fam, .log extensions)
ENSG00000116704.7_blood (with .bed, .bim, .fam, .log extensions)
.bed is the genotype file (you cannot view it unless you load it into R)
.bim is the list of SNPs in the cis window for that gene (+/- 500kb of the TSS)
.fam is the list of people for that gene/tissue pair which records the gene expression values in column 6 and the people IDs in column 2 
.log is the log file from when I created these files 

Our first exercise is to estimate heritability of the two genes above in each of the two tissues listed above. 
You will use FUSION to do this, which runs GCTA (the standard method for estimating heritability). (http://gusevlab.org/projects/fusion/)

In command line: (you should be in the BGGN239/ directory)

```
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
rm master.zip
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip
rm master.zip
```

If on Expanse, launch an interactive session and load R: 
```
module load gcc/9.2.0
module load r
```
### Launch R and install required libraries:
```
R
> install.packages(c('optparse','RColorBrewer'))
> install.packages("Rcpp")
> install.packages("Matrix")
> install.packages("RcppEigen")
> install.packages('plink2R-master/plink2R/',repos=NULL)
#> install.packages("glmnet") #if problem with install, simply comment out library(glmnet) in fusion_twas-master/FUSION.compute_weights.R
```

### Download plink
```
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip 
#rename LICENSE to LICENSE_plink when prompted 
./plink #should run 
```

For example, run this in command line:

```
gene=ENSG00000158864.12
tissue=brain
mkdir tmp
home=XX #PUT THE FULL PATH UP TO BGGN239 REPO (see for PATH_plink)
Rscript fusion_twas-master/FUSION.compute_weights.R --bfile ${gene}_${tissue} --covar covar_${tissue}.txt --tmp tmp/tmp_${gene} --out out_${gene}_${tissue} --models top1 --PATH_gcta fusion_twas-master/gcta_nr_robust --PATH_plink ${home}/BGGN239/plink --hsq_p 1 >> gcta_${gene}_${tissue}.txt
````

Running this script also produces a file: out_ENSG00000158864.12_brain.wgt.RDat. We will use this file later. 

For the first in-class exercise, determine the GCTA-estimated heritability of each gene (n = 2) in each tissue (n = 2), the heritability standard error, and the heritability p-value. 
Which genes in which tissues have significantly positive heritability?
Are the heritability estimates significantly different from one another? What is the p-value for the difference? 
 
For the second in-class exercise, use the read_plink() R function to read in genotypes for ENSG00000158864.12 in blood. You can look up an example of this function as it's used in fusion_twas-master/FUSION.compute_weights.R. Then, calculate the squared correlation of genotypes across all pairs of SNPs. How many pairs of variants have r2 = 1, r2 >= 0.8? What proportion of SNP pairs are in perfect LD (e.g. r2 = 1). 

For the third in-class exercise, we will perform a TWAS: 

```
for sumstats in PASS_Alzheimers_deRojas2021.sumstats PASS_UC_deLange2017.sumstats
do 
Rscript fusion_twas-master/FUSION.assoc_test.R --sumstats $sumstats --weights twas.pos --weights_dir ${home}/BGGN239 --ref_ld_chr 1000G.EUR. --chr 1 --out TWAS_${sumstats}.dat
done
```

What are the TWAS p-values for each gene in each tissue? Which genes are significantly associated with each disease? 


