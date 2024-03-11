## Overview
ParaMask encompasses three main steps:

### 1. PrepareVCF_fromVCF
In this initial phase, a VCF (Variant Call Format) file undergoes processing to yield three ParaMaskInput files:

- **.het** Contains essential data for the Expectation Maximization algorithm and Read Ratio Deviations.
- **.cov.stat.txt** Encodes coverage per sample and per site.
- **.cov.gw.txt** Displays genome-wide coverage of non-missing sites per individual.

### 2. ParaMask_EM
This stage harnesses population genomic signatures of multicopies for seed SNP generation, employing the following procedures:

- **EM algorithm:** Simultaneously fits two Beta-binomial regressions on heterozygote frequency as a function of the Minor Allele frequencies. One regression pertains to single-copy regions, while the other addresses multicopy regions. Classification hinges on the Log-Likelihood Ratio (LLR).
- **RRD testing:** Utilizes the mean and variance from the read ratio deviation of single-copy classified SNPs to construct a normal confidence interval, validating previously uncertain classified SNPs.
- **EM algorithm for distance dissection:** Pinpoints distances between seed SNPs within and between multicopy regions. It calculates mean parameters using a mixture of geometric distributions. The cutoff distance is established where the two geometrics have identical density. For increased stability, this process defaults to 1000 repetitions, with the median cutoff selected.

Creates two intermediate Output files

- **.EMresults.het** An updated het file with het with stats on EM classification
- **.dist** File containing the distance cutoff

Automated plotting generates visualization:



### 3. ParaMask_Cluster_Seeds
In the final step, SNPs are clustered into multicopy haplotypes, and a comprehensive SNP annotation is provided. This phase yields three output files:

- **.finalClass.het** file The original Het file with the definitive status.
- **.clusters.txt** Cluster file outlining each multicopy SNP along with its annotation.
- **.finalClass.bed** Bed file distinguishing single- and multicopy regions.

## installation

### Download
- ```bash
  git clone https://github.com/Fulgione-group/ParaMask.git

### Prerequisite

**tested on**
#### R
- R version 4.0.4
- VGAM version 1.1.8
- ggplot2 version 3.4.0

#### Java
- byte code generated for java version   1.8 using gradle
- You can compile bytecode from source located in the ParaMask/src folder

## Usage 

### 1. PrepareVCF_fromVCF
- ```bash
  #!/bin/bash
  java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/PrepareParaMaskInput_fromVCF.jar\
        --vcf/-v $INPUT_VCF #Input VCF
  # Optional parameters
        --missingness/-m $MAX_MISSINGNESS_PROPORTION  #float, default = 0: no missing sites allowed
        --popfile/-p $PATH_TO_POPFILE #Full path to popfile, a list of samples in each row, default all samples in the VCF
        --out/-o $PATH_TO OUTPUT #full path to the Output file, dafault is the input file. Extensions for the different files are added automatically
        --noVaryingFormat/nVF #NOT RECOMMENDED sets Varying genotype format of the VCF to false, default true.

### 2. ParaMask_EM
- ```bash
  #!/bin/bash
  Rscript --vanilla $PATH_TO_INSTALLATION_FOLDER/ParaMask/ParaMask_EM_v2.4.R\
        --hetfile/-h $PATH_TO_HET_FILE
  # Optional parameters
        --outpath/-o $PATH_TO_OUTPUTDIR
        --missingness/-m $MAX_MISSINGNESS_PROPORTION  #float, default = 0.1: no missing sites allowed
        --verbose/-v #Verbose shows current steps of ParaMask, fitting process of VGAM, default is false
        --ID $RUN_ID #ID for ParaMask_EM_run, will be used in file naming
        --chrom/-c $CHROMOSOME #Use a specific chromosome only, default: all chromosomes
        --noRRD #do not use RRD to classify uncertain, defaut: True
        --tolerance/-t $EM_TOLERANCE$ #Tolerance for Parameters estimated by the EM algorithm on heterezygote frequency, default: 0.001
        --startline/-s $INTEGER #If you want to analyses a certain subset of SNPs in the hetfile you can specify start end lines
        --endline/-e $INTEGER #If you want to analyses a certain subset of SNPs in the hetfile you can specify start end lines
        --boundary/-b $FLOAT #NOT RECOMMENDED effectively constraints the upper Parameter space of the MAF*(Z=="K") variable, Helps with EM convergence in extreme cases


### 3. ParaMask_Cluster_Seeds
- ```bash
   #!/bin/bash
   java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/ParaMask_Cluster_Seeds.jar\
        --cov PATH_TO_COV_FILE/-c $PATH_TO_COVSTAT_FILE #
        --het PATH_TO_HET_FILE/-h $PATH_TO_HET_FILE
        --covgw PATH_TO_COVGW_FILE/-cg $PATH_TO_COVGW_FILE
        --cutoff/-cd $INTEGER #distance cutoff
        --range $INTEGER,INTEGER #CHR_START,CHR_END
        --purge $INTEGER #Clusters with number of SNPs <=INTEGER are purged, default = 1.


## Output files

### .finalClass.het

This file contains all per SNP statistics and results used for classification

Columns:
1. Chromosome
2. Position
3. Number of non missing genotypes
4. Minor allele frequency
5. Heterozygous genotype frequency
6. Homozygous genotype frequency of the major allele
7. Homozyhous genotype frequency of the minor allele
8. Mean coverage across genotypes
9. Mean coverage across homozygous genotypes
10. Mean coverage across heterozygote genotypes
11. Reference allele depth of all heterozygote genotypes
12. Alternative allele depth
13. Allelic ratio (11/12)
14. Read ratio deviation (RRD, based on 13)
15. Likelihood of SNP beeing single-copy
16. Likelihood of SNP beeing multicopy
17. Likelihood ratio
18. Classification after EM step: 0 = single-copy; 1 = uncertain; 2 = multicopy-Seed
19. Seed because of allelic ratio deviation: 0 = no Seed/ Seed based on EM; 1 = Seed because of RRD (EM classified uncertain before)
20. Final status: 0= single copy; 1 = muliticopy
21. cluster Number: 0 = no cluster (single copy); 1...N = multicopy

| Chromosome | Position | Non.missing | Minor.allele.freq | Heterozygous.geno.freq | Homozygous1.geno.freq | Homozygous2.geno.freq | Mean.coverage | Mean.coverage.hom | Mean.coverage.het | Het.reference.allele.depth | Het.alt.allele.depth | Het.allele.ratio | Het.allele.deviation | L1 | L2 | LLR | EM_class | allele.deviation.seed | finalClass | cluster | 
|------------|----------|-------------|---------------------|-------------------------|-------------------------|-------------------------|----------------|---------------------|---------------------|---------------------------|-----------------------|-------------------|-----------------------|----|----|-----|----------|------------------------|------------|---------|
| chr1       | 11065    | 100         | 0.025               | 0.05                    | 0.95                    | 0                       | 19.96          | 19.83158             | 22.4                | 80                        | 32                    | 0.71428573        | 4.5355735 | 0.411559811962468 |  0.963941228136795 | -0.851075965471138	      | 2  | 1  |  1          | 1       |
| chr1       | 11226    | 99          | 0.26767677          | 0.35353535              | 0.5555556               | 0.09090909                       | 10.686869      | 10.453125            | 11.114285           | 193                       | 196                   | 0.49614397        | -0.15210603 | 0.0904636982907597 |  7.99716717457639e-10  | 18.5439569215409          | 0  | 0  | 0          | 0       |
| chr1       | 11612    | 98          | 0.4489796           | 0.48979592              | 0.30612245              | 0.2040816                     | 11.397959      | 10.5                | 12.333333           | 312                       | 280                   | 0.527027           | 1.3151919 |  0.0562953395595597 | 8.33420025699269e-16 | 31.8438504003201            | 0  | 0  | 0          | 0       |
| chr1       | 11993    | 99          | 0.43434343         | 0.4040404               | 0.36363637              | 0.2040816                      | 9.474748       | 9.338983             | 9.675               | 185                       | 202                   | 0.47803617        | -0.8641586 | 0.0423081470021819	| 6.12813943557946e-19	| 38.773449969909             | 0  | 0  | 0          | 0       |
| chr1       | 12373    | 99          | 0.030303031         | 0.060606062             | 0.93939394              | 0                       | 10.565657      | 10.634409            | 9.5                 | 30                        | 27                    | 0.5263158         | 0.3973597 | 0.338789194072882 | 0.957120910950128 | -1.03855165964169 | 1  | 0  | 0          | 0       |




### .clusters.txt

This files contains additional per multicopy SNP statistics, with details on why they classified and which genotypes are involved.

Columns:
1. Chromosome
2. Position
3. Cluster Number
4. Classification after EM step: 0 = single-copy; 1 = uncertain; 2 = multicopy-Seed
5. Reason for multicopy SNP classification: seed = Seed SNP; hetcov = single-copy SNP in multicopy region classified by coverage (at current position) of heterozygote genotypes at the last Seed; bridge = uncertain SNP with no excess of coverage, but in between 2 SNPs with multicopy signals
6. Coverage (at current position) of heterozygote genotypes at the last Seed
7. Heterozygote genotypes at the last Seed: Colon seperated list. Can be used to extract genotype specific Clusters

| Chromosome | Position | Cluster | EmClass | ClusterCause | CovHetOfLastSeed | CovGWHetOfLastSeed | HetGenOfLastSeed                                               |
|------------|----------|---------|---------|--------------|------------------|---------------------|---------------------------------------------------------------|
| chr1       | 10491    | 1       | 2       | seed         | 18.465643        | 12.807346           | genotype_1:genotype_2:genotype_3:genotype_4:genotype_5:genotype_6:genotype_7:genotype_8:genotype_9:genotype_10:genotype_11:genotype_13:genotype_14:genotype_17:genotype_18:genotype_19:genotype_20:genotype_22:genotype_23:genotype_24:genotype_25:genotype_28:genotype_30:genotype_32:genotype_33:genotype_35:genotype_37:genotype_38:genotype_39:genotype_40:genotype_42:genotype_43:genotype_44:genotype_45:genotype_46:genotype_48:genotype_50:genotype_51:genotype_52:genotype_53:genotype_54:genotype_55:genotype_56:genotype_57:genotype_58:genotype_59:genotype_60:genotype_61:genotype_63:genotype_64:genotype_65:genotype_67:genotype_68:genotype_69:genotype_70:genotype_71:genotype_72:genotype_73:genotype_74:genotype_75:genotype_76:genotype_77:genotype_78:genotype_80:genotype_81:genotype_82:genotype_84:genotype_86:genotype_92:genotype_93:genotype_94:genotype_95:genotype_96:genotype_97:genotype_98:genotype_99 |
| chr1       | 10492    | 1       | 2       | seed         | 21.133333        | 12.640491          | genotype_12:genotype_54:genotype_55:genotype_56:genotype_65:genotype_72:genotype_77:genotype_88 |
| chr1       | 10569    | 1       | 2       | seed         | 18.125           | 12.657732          | genotype_2:genotype_3:genotype_4:genotype_5:genotype_6:genotype_7:genotype_8:genotype_9:genotype_10:genotype_11:genotype_13:genotype_14:genotype_15:genotype_16:genotype_17:genotype_18:genotype_19:genotype_20:genotype_22:genotype_23:genotype_24:genotype_25:genotype_28:genotype_30:genotype_31:genotype_32:genotype_33:genotype_35:genotype_37:genotype_38:genotype_39:genotype_40:genotype_42:genotype_43:genotype_44:genotype_48:genotype_49:genotype_50:genotype_51:genotype_52:genotype_53:genotype_54:genotype_55:genotype_56:genotype_57:genotype_58:genotype_59:genotype_60:genotype_61:genotype_63:genotype_64:genotype_65:genotype_68:genotype_69:genotype_70:genotype_71:genotype_73:genotype_74:genotype_75:genotype_76:genotype_77:genotype_78:genotype_79:genotype_80:genotype_81:genotype_82:genotype_83:genotype_84:genotype_86:genotype_87:genotype_88:genotype_89:genotype_92:genotype_93:genotype_94:genotype_95:genotype_96:genotype_97:genotype_98:genotype_99 |
| chr1       | 10593    | 1       | 1       | hetcov       | 19.7125          | 12.643446          | genotype_2:genotype_3:genotype_4:genotype_5:genotype_6:genotype_7:genotype_8:genotype_9:genotype_10:genotype_11:genotype_13:genotype_14:genotype_15:genotype_16:genotype_17:genotype_18:genotype_19:genotype_20:genotype_22:genotype_23:genotype_24:genotype_25:genotype_28:genotype_30:genotype_31:genotype_32:genotype_33:genotype_35:genotype_37:genotype_38:genotype_39:genotype_40:genotype_42:genotype_43:genotype_44:genotype_48:genotype_49:genotype_50:genotype_51:genotype_52:genotype_53:genotype_54:genotype_55:genotype_56:genotype_57:genotype_58:genotype_59:genotype_60:genotype_61:genotype_63:genotype_64:genotype_65:genotype_68:genotype_69:genotype_70:genotype_71:genotype_73:genotype_74:genotype_75:genotype_76:genotype_77:genotype_78:genotype_79:genotype_80:genotype_81:genotype_82:genotype_83:genotype_84:genotype_86:genotype_87:genotype_88:genotype_89:genotype_92:genotype_93:genotype_94:genotype_95:genotype_96:genotype_97:genotype_98:genotype_99 |
| chr1       | 10641    | 1       | 2       | seed         | 19.525           | 12.643446          | genotype_5:genotype_15:genotype_68:genotype_78:genotype_93       |


### .finalClass.bed

**For most users probably the most important output file**
1-based bed bed file containing genomic regions and copy number status (single / multicopy), number of SNPs and Cluster Number.

Columns:
1. Chromosome
2. Start of genomic region
3. End of the genomic region
4. Region status: 0 = single-copy; 1 = multicopy
5. Number of SNPs
6. Cluster Number

| Chromosome | Start | End   | Type:0-single-copy; 1-multi-copy) | nSNPs | Cluster |
|------------|-------|-------|----------------------------------------|------|---------|
| chr1       | 1     | 10473 | 0                                      | 63   | 0       |
| chr1       | 10474 | 11145 | 1                                      | 12   | 1       |
| chr1       | 11146 | 19264 | 0                                      | 45   | 0       |
| chr1       | 19265 | 22208 | 1                                      | 51   | 2       |
| chr1       | 22209 | 36211 | 0                                      | 80   | 0       |


## Example files


### Run forward simulations with SeDuS


### Process Output files
