\***Attention***

## Change log

### 04.September.2025

Bugs with scientific formating and the missingness filter have been corrected in the ParaMaskEM package. A new test version 1.0.1 has been uploaded with a new flag "--nSNPs [int]" for fitting the EM with a reduced number of [int] randomly selected SNPs, while providing classification of all SNPs.

### 24.April.2025

The EM step has been modularized into an R package to improve code organization and streamline dependency management. The package is located in the ParaMaskEM/ directory, and installation instructions are provided [here](#install-r-package-paramaskem-from-within-r).

### June.2014
A new EM intialization has been implemented in version ParaMask_EM_v0.2.7.1.R, version 0.2.5 is deprecated!!
The latest version ParaMask_EM_v0.2.7.2.R fixes some issues in the LLR calculation for low sample sizes.


## Overview
<br>

ParaMask is a method to identify multicopy regions in population-level whole genome data, including tandem or segmental duplications, copy number variants, gene families with various paralog copies, transposable elements, and other repeats. 
To run ParaMask, we first prepare input files from a vcf (script 1), then we run a first classification of SNPs based on excess heterozygosity and read ratio deviations (script 2), then we cluster collapsed SNPs in multicopy haplotypes (script 3).
For questions please contact us via [email](btjeng@mpipz.mpg.de). 

### 1. PrepareParaMaskInput_fromVCF
\*Note: please do not filter for excess of coverage, since this signal is used in the ParaMask haplotype clustering algorithm.\
First, we process a VCF (Variant Call Format) file containing SNPs and output three ParaMaskInput files:


- **.het** Contains essential data for the Expectation Maximization algorithm and Read Ratio Deviations.
- **.cov.stat.txt** Encodes coverage per sample and per site.
- **.cov.gw.txt** Displays genome-wide coverage of non-missing sites per individual.

<br>

### 2. ParaMask_EM
This step produces a first classification of SNPs in single-copy and multicopy regions, employing the following procedures:

- **EM algorithm:** Simultaneously fits two Beta-binomial regressions on heterozygote frequency as a function of the minor allele frequencies. One regression fits single-copy regions, while the other fits multicopy regions. SNPs are classified based on the Log-Likelihood Ratio (LLR).
- **RRD testing:** Utilizes the mean and variance from the read ratio deviation of SNPs classifies as single-copy in the EM-step, to construct a normal confidence interval. This improves the power to detect SNPs that were classified as uncertain in the EM-step.
- **EM algorithm for distance dissection:** This step fits a mixture of two geometric distributions to the distances among seed SNPs, one for distances within regions and one for distances between regions. The cutoff distance is established where the two geometrics have identical density. For increased stability, by default this process is repeated 1000 times, and the distance cutoff is set to the median across replicates.

This step automatically generates plots for visualization (see diagnostic plots in the examples section), and it creates two intermediate output files:

- **.EMresults.het** An updated het file with statistics on the EM classification (detail below)
- **.dist** File containing the distance cutoff


<br>

### 3. ParaMask_Cluster_Seeds
In the final step, SNPs are clustered into multicopy haplotypes, and SNPs are classified into single- and multicopy SNPs. If multiple chromosomes are present, this step needs to be run seperately on them. This step outputs the following files:

- **.finalClass.het** The original Het file with the final classification of SNPs.
- **.clusters.txt** Cluster file outlining each multicopy SNP along with its annotation.
- **.finalClass.bed** Bed file with single- and multicopy regions.

<br>
<br>
<br>
<br>
<br>
<br>

## installation

<br>

### Download
- ```bash
  git clone https://github.com/Fulgione-group/ParaMask.git

### Install R package ParaMaskEM from within R
```
# Install devtools if needed
install.packages("devtools")

# Install the ParaMaskEM package from your repo subdirectory
devtools::install_github("Fulgione-group/ParaMask", subdir = "ParaMaskEM")

# locate script to run the pipeline with command line args
system.file("scripts", "run_ParaMask_EM.R", package = "ParaMaskEM")
```

### Prerequisite

<br>

#### R
- R version >= 4.0.4
- VGAM version >= 1.1.1
- ggplot2 version >= 3.4.0
- patchwork >= 1.0.0
- data.table >= 1.12.8

#### Java
- java version 1.8 

<br>
<br>
<br>
<br>
<br>
<br>

## Usage 

<br>

### 1. PrepareVCF_fromVCF

Example:
```bash
#!/bin/bash
java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/PrepareParaMaskInput_fromVCF.jar\
        --vcf $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf\
        --missingness 0.1

```
<br>

| Option                          | Description                                                |
|---------------------------------|------------------------------------------------------------|
|**Required**|
| **--vcf/-v**                    | Input VCF                                      |
|**Optional**|
| **--missingness/-m**            | Maximum proportion of missing genotypes per site. Default = 0: no missing sites allowed |
| **--popfile/-p**                | Full path to popfile, a list of samples in each row. Default NULL, all samples in the VCF |
| **--out/-o**                    | Full path to the output file, default is the input file. Extensions for the different files are added automatically |


<br>

### 2. ParaMask_EM

Example:
```bash
#!/bin/bash
Rscript --vanilla $PATH_to_pipeline_script_from_ParaMaskEM/run_ParaMask_EM.R\
        --het $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.het.stat.txt\
        --missingness 0.1\
        --outdir $PATH_TO_INSTALLATION_FOLDER/Example_files/Output\
        --ID test  
```
<br>
   
| Option                | Description |
|-----------------------|-------------|
|**Required**|
| **--het/-h**      | Input full path to het file |
|**Optional**|
| **--outdir/-o**      | Input full path to the output directory |
| **--missingness/-m**  | Input float, default = 0: no missing sites allowed |
| **--verbose/-v**      | Verbose shows current steps of ParaMask, fitting process of VGAM, default is false |
| **--ID**      | Input ID for file naming |
| **--chrom/-c**        | Input chromosome name to only use a specific chromosome. Default: all chromosomes |
| **--noRRD**           | Do not use read ratio deviations. Default: True |
| **--tolerance/-t**    | Input tolerance for parameters estimated by the EM algorithm on heterozygote frequency, default: 0.001 |
| **--startline/-s**    | Integer: Starting line of the het file. Default=2 |
| **--endline/-e**      | Integer: Ending line of the het file. Default last line|
| **--boundary/-b**     | **NOT RECOMMENDED** Float: constrain to the lower,upper limit of the MAF*(Z=="K") parameter. This can help with EM convergence in cases where SNPs are clustered in a small range of maf. If boundaries are exceeded a modified step takes with the lower or upper limit as offset is taken. Disabled by default.|


<br>

### 3. ParaMask_Cluster_Seeds

Example:
```bash

#!/bin/bash
java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/ParaMask_Cluster_Seeds.jar\
        --cov $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.cov.stat.txt\
        --het $PATH_TO_INSTALLATION_FOLDER/Example_files/Output/test_EMresults.het\
        --covgw  $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.cov.gw.txt\
        --cutoff  $(tail -1 $PATH_TO_INSTALLATION_FOLDER/Example_files/Output/test_EMresults.dist)\
        --range 1,1000000
```

<br>

| Option                                  | Description                                                  |
|-----------------------------------------|--------------------------------------------------------------|
| **Required**|
| **--cov/-c**           | Full path to coverage file per sample and per site |
| **--het/-h**           | Full path to het file generated in step 2 |
| **--covgw/-cg**      | Full path to genome-wide coverage of non-missing sites per individual |
| **--cutoff/-cd**                        | Integer: Distance cutoff                                     |
| **--range**                             | Integer,Integer: CHR_START,CHR_END                           |
|**Optional**|
| **--purge**                             | Integer: Multicopy regions with at most this specified number of SNPs are purged. Default = 1 |

<br>
<br>
<br>
<br>
<br>
<br>


## Output files

<br>

### .finalClass.bed

**For most users probably the most important output file** <br>
1-based bed file containing genomic regions and copy number status (single- / multicopy), number of SNPs and cluster number.

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

<br>

### .finalClass.het


This file contains all statistics per SNPs and the results used for classification

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
12. Alternative allele depth of all heterozygote genotypes
13. Allelic ratio (11/12)
14. Read ratio deviation (RRD, based on 13)
15. Likelihood of SNP beeing single-copy
16. Likelihood of SNP beeing multicopy
17. Log-Likelihood-Ratio
18. Classification after EM step: 0 = single-copy; 1 = uncertain; 2 = multicopy-Seed
19. Reason why a SNP is classified as seed: 0 = no Seed/ Seed based on EM; 1 = Seed because of RRD (EM classified uncertain before)
20. Final status: 0= single copy; 1 = muliticopy
21. cluster Number: 0 = no cluster (single copy); 1...N = multicopy

| Chromosome | Position | Non.missing | Minor.allele.freq | Heterozygous.geno.freq | Homozygous1.geno.freq | Homozygous2.geno.freq | Mean.coverage | Mean.coverage.hom | Mean.coverage.het | Het.reference.allele.depth | Het.alt.allele.depth | Het.allele.ratio | Het.allele.deviation | L1 | L2 | LLR | EM_class | allele.deviation.seed | finalClass | cluster | 
|------------|----------|-------------|---------------------|-------------------------|-------------------------|-------------------------|----------------|---------------------|---------------------|---------------------------|-----------------------|-------------------|-----------------------|----|----|-----|----------|------------------------|------------|---------|
| chr1       | 11065    | 100         | 0.025               | 0.05                    | 0.95                    | 0                       | 19.96          | 19.83158             | 22.4                | 80                        | 32                    | 0.71428573        | 4.5355735 | 0.411559811962468 |  0.963941228136795 | -0.851075965471138	      | 2  | 1  |  1          | 1       |
| chr1       | 11226    | 99          | 0.26767677          | 0.35353535              | 0.5555556               | 0.09090909                       | 10.686869      | 10.453125            | 11.114285           | 193                       | 196                   | 0.49614397        | -0.15210603 | 0.0904636982907597 |  7.99716717457639e-10  | 18.5439569215409          | 0  | 0  | 0          | 0       |
| chr1       | 11612    | 98          | 0.4489796           | 0.48979592              | 0.30612245              | 0.2040816                     | 11.397959      | 10.5                | 12.333333           | 312                       | 280                   | 0.527027           | 1.3151919 |  0.0562953395595597 | 8.33420025699269e-16 | 31.8438504003201            | 0  | 0  | 0          | 0       |
| chr1       | 11993    | 99          | 0.43434343         | 0.4040404               | 0.36363637              | 0.2040816                      | 9.474748       | 9.338983             | 9.675               | 185                       | 202                   | 0.47803617        | -0.8641586 | 0.0423081470021819	| 6.12813943557946e-19	| 38.773449969909             | 0  | 0  | 0          | 0       |
| chr1       | 12373    | 99          | 0.030303031         | 0.060606062             | 0.93939394              | 0                       | 10.565657      | 10.634409            | 9.5                 | 30                        | 27                    | 0.5263158         | 0.3973597 | 0.338789194072882 | 0.957120910950128 | -1.03855165964169 | 1  | 0  | 0          | 0       |


<br>

### .clusters.txt



This files contains additional statistics per multicopy SNPs, with details on why they were classified as multicopy and which genotypes are involved.

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

<br>

### Diagnostic plots

ParaMask outputs diagnostic plots in pdf format

1. **iterationN.pdf**
  * diagnostic plots of posterioir weights of the EM after the Nth iteration
2. **LLR.pdf**
  * Plot of the Log-Likelihood-Ratio
3. **AR.pdf** and **RRD.pdf**
  * Plots of allelic ratios and read ratio deviations (and densities) grouped by EM classfication (single-copy, multicopy, uncertain)
4. **dist.pdf**
  * Distance plot
<br>
<br>
<br>
<br>
<br>
<br>

## Example files

Can be reproduced by running example scripts above, also located in the scripts folder

Folders:
1. **Input** 
* **Simulations_10PercentDuplications.vcf** Simulated VCF with 10% duplications
* **Simulated_Regions.bed** Simulated regions
* **stat.txt**, **cov.gw.txt** and **cov.stat.txt** Computed files by PrepareParaMaskInput_fromVCF
2. **Output**
  * All files generated by Step 2 and 3.

