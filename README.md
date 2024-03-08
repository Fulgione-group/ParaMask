# Overview
ParaMask encompasses three main steps:

## 1. PrepareVCF_fromVCF
In this initial phase, a VCF (Variant Call Format) file undergoes processing to yield three ParaMaskInput files:

- **.het** Contains essential data for the Expectation Maximization algorithm and Read Ratio Deviations.
- **.cov.stat.txt** Illustrates coverage per sample and per site.
- **.cov.gw.txt** Displays genome-wide coverage of non-missing sites per individual.

## 2. ParaMask_EM
This stage harnesses population genomic signatures of multicopies for seed SNP generation, employing the following procedures:

- **EM algorithm:** Simultaneously fits two Beta-binomial regressions on heterozygote frequency as a function of the Minor Allele frequencies. One regression pertains to single-copy regions, while the other addresses multicopy regions. Classification hinges on the Log-Likelihood Ratio (LLR).
- **RRD testing:** Utilizes the mean and variance from the read ratio deviation of single-copy classified SNPs to construct a normal confidence interval, validating previously uncertain classified SNPs.
- **EM algorithm for distance dissection:** Pinpoints distances between seed SNPs within and between multicopy regions. It calculates mean parameters using a mixture of geometric distributions. The cutoff distance is established where the two geometrics have identical density. For increased stability, this process defaults to 1000 repetitions, with the median cutoff selected.

## 3. ParaMask_Cluster_Seeds
In the final step, SNPs are clustered into multicopy haplotypes, and a comprehensive SNP annotation is provided. This phase yields three output files:

- **.finalClass.het** file The original Het file with the definitive status.
- **.clusters.txt** Cluster file outlining each multicopy SNP along with its annotation.
- **.finalClass.bed** Bed file distinguishing single- and multicopy regions.

# Details 

## 1. PrepareVCF_fromVCF
- ```bash
  #!/bin/bash
  java -jar /PATH_TO_PARAMASK_BYTECODE/PrepareParaMaskInput_fromVCF.jar\
        --vcf $INPUT_VCF #Input VCF
  # Optional parameters
        --missingness $MAX_MISSINGNESS_PROPORTION  #float, default = 0: no missing sites allowed
        --popfile $PATH_TO_POPFILE #Full path to popfile, a list of samples in each row, default all samples in the VCF
        --out $PATH_TO OUTPUT #full path to the Output file, dafault is the input file. Extensions for the different files are added automatically
        --noVaryingFormat #Not recommended sets Varying genotype format of the VCF to false, default true.

## 1. ParaMask_EM
- ```bash
  Rscript --vanilla ~/PATH_TO_PARAMASK_SRC/ParaMask_EM_v2.4.R\
        --hetfile $PATH_TO_HET_FILE
  # Optional parameters
        --outpath $PATH_TO_OUTPUTDIR
        --missingness $MAX_MISSINGNESS_PROPORTION  #float, default = 0.1: no missing sites allowed
        --verbose #Verbose shows current steps of ParaMask, fitting process of VGAM, default is false
        --ID $RUN_ID #ID for ParaMask_EM_run, will be used in file naming
        --chrom $CHROMOSOME #Use a specific chromosome only, default: all chromosomes
        --noRRD #do not use RRD to classify uncertain, defaut: True
  




### Run forward simulations with SeDuS


### Process Output files
