## Overview
ParaMask encompasses three main steps:

### 1. PrepareVCF_fromVCF
In this initial phase, a VCF (Variant Call Format) file undergoes processing to yield three ParaMaskInput files:

- **.het** Contains essential data for the Expectation Maximization algorithm and Read Ratio Deviations.
- **.cov.stat.txt** Illustrates coverage per sample and per site.
- **.cov.gw.txt** Displays genome-wide coverage of non-missing sites per individual.

### 2. ParaMask_EM
This stage harnesses population genomic signatures of multicopies for seed SNP generation, employing the following procedures:

- **EM algorithm:** Simultaneously fits two Beta-binomial regressions on heterozygote frequency as a function of the Minor Allele frequencies. One regression pertains to single-copy regions, while the other addresses multicopy regions. Classification hinges on the Log-Likelihood Ratio (LLR).
- **RRD testing:** Utilizes the mean and variance from the read ratio deviation of single-copy classified SNPs to construct a normal confidence interval, validating previously uncertain classified SNPs.
- **EM algorithm for distance dissection:** Pinpoints distances between seed SNPs within and between multicopy regions. It calculates mean parameters using a mixture of geometric distributions. The cutoff distance is established where the two geometrics have identical density. For increased stability, this process defaults to 1000 repetitions, with the median cutoff selected.

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



## Details 

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
- ```
  Chromosome	Position	Non.missing	Minor.allele.freq	Heterozygous.geno.freq	Homozygous1.geno.freq	Homozygous2.geno.freq	Mean.coverage	Mean.coverage.hom	Mean.coverage.het	Het.reference.allele.depth	Het.alt.allele.depth	Het.allele.ratio  Het.allele.deviation	L1	L2	LLR	EM_class	allele.deviation.seed	finalClass	cluster
  chr1	10641	100	0.025	0.05	0.95	0	19.68	19.68421	19.6	72	26	0.7346939	4.646702	0.411559811962468	0.963941228136795	-0.851075965471138	2	1 1	1
  chr1	10680	100	0.04	0.08	0.92	0	19.63	19.554348	20.5	118	46	0.7195122	5.6222553	0.22685182553148	0.943842011392948	-1.42566173794247	2	1 1	1
  chr1	10776	99	0.08080808	0.14141414	0.8484849	1	20.818182	20.964706	19.928572	185	94	0.6630824	5.4480276	0.203133317411145	0.0113927604667986	2.88088439259977	2	1	1	1
  chr1	10798	100	0.49	0.94	0.02	4	20.9	14	21.340425	993	1013	0.49501497	-0.4465443	2.90395537401198e-14	0.0210172995916636	-27.3077081818396	2	0  1	1
  chr1	10811	99	0.05050505	0.1010101	0.8989899	0	19.959597	19.853933	20.9	143	66	0.68421054	5.3262014	0.147898455696237	0.931023048035841	-1.83975810506303  2	1  1	1
  chr1	10906	100	0.015	0.03	0.97	0	19.87	19.927835	18	38	16	0.7037037	2.993821	0.596228304224322	0.977961587648625	-0.494846738335215	2  1  1  1
  chr1	10971	99	0.5	1	0	0	21.424242	0	21.424242	1063	1058	0.5011787	0.10856746	2.3785162706989e-19	0.608498643720493	-42.3858792938336	2  0	1  1
  chr1	11065	100	0.025	0.05	0.95	0	19.96	19.83158	22.4	80	32	0.71428573	4.5355735	0.411559811962468	0.963941228136795	-0.851075965471138	2  1  1  1
  chr1	11226	99	0.26767677	0.35353535	0.5555556	9	10.686869	10.453125	11.114285	193	196	0.49614397	-0.15210603	0.0904636982907597	7.99716717457639e-10	18.5439569215409  0  0  0  0
  chr1	11612	98	0.4489796	0.48979592	0.30612245	20	11.397959	10.5	12.333333	312	280	0.527027	1.3151919	0.0562953395595597	8.33420025699269e-16	31.8438504003201	0  0  0  0

### .finalClass.het
```markdown
 Chromosome | Position | Non.missing | Minor.allele.freq | Heterozygous.geno.freq | Homozygous1.geno.freq | Homozygous2.geno.freq | Mean.coverage | Mean.coverage.hom | Mean.coverage.het | Het.reference.allele.depth | Het.alt.allele.depth | 
 Het.allele.ratio | Het.allele.deviation | L1 | L2 | LLR | EM_class | allele.deviation.seed | finalClass | cluster
 --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
 chr1 | 10641 | 100 | 0.025 | 0.05 | 0.95 | 0 | 19.68 | 19.68421 | 19.6 | 72 | 26 | 0.7346939 | 4.646702 | 0.411559811962468 | 0.963941228136795 | -0.851075965471138 | 2 | 1 1 | 1
 chr1 | 10680 | 100 | 0.04 | 0.08 | 0.92 | 0 | 19.63 | 19.554348 | 20.5 | 118 | 46 | 0.7195122 | 5.6222553 | 0.22685182553148 | 0.943842011392948 | -1.42566173794247 | 2 | 1 1 | 1
 chr1 | 10776 | 99 | 0.08080808 | 0.14141414 | 0.8484849 | 1 | 20.818182 | 20.964706 | 19.928572 | 185 | 94 | 0.6630824 | 5.4480276 | 0.203133317411145 | 0.0113927604667986 | 2.88088439259977 | 2 | 1 | 1 | 1
 chr1 | 10798 | 100 | 0.49 | 0.94 | 0.02 | 4 | 20.9 | 14 | 21.340425 | 993 | 1013 | 0.49501497 | -0.4465443 | 2.90395537401198e-14 | 0.0210172995916636 | -27.3077081818396 | 2 | 0 1 | 1
 chr1 | 10811 | 99 | 0.05050505 | 0.1010101 | 0.8989899 | 0 | 19.959597 | 19.853933 | 20.9 | 143 | 66 | 0.68421054 | 5.3262014 | 0.147898455696237 | 0.931023048035841 | -1.83975810506303 | 2 | 1 1 | 1
 chr1 | 10906 | 100 | 0.015 | 0.03 | 0.97 | 0 | 19.87 | 19.927835 | 18 | 38 | 16 | 0.7037037 | 2.993821 | 0.596228304224322 | 0.977961587648625 | -0.494846738335215 | 2 | 1 1 1 | 1
 chr1 | 10971 | 99 | 0.5 | 1 | 0 | 0 | 21.424242 | 0 | 21.424242 | 1063 | 1058 | 0.5011787 | 0.10856746 | 2.3785162706989e-19 | 0.608498643720493 | -42.3858792938336 | 2 | 0 | 1 1
 chr1 | 11065 | 100 | 0.025 | 0.05 | 0.95 | 0 | 19.96 | 19.83158 | 22.4 | 80 | 32 | 0.71428573 | 4.5355735 | 0.411559811962468 | 0.963941228136795 | -0.851075965471138 | 2 | 1 1 1 | 1
 chr1 | 11226 | 99 | 0.26767677 | 0.35353535 | 0.5555556 | 9 | 10.686869 | 10.453125 | 11.114285 | 193 | 196 | 0.49614397 | -0.15210603 | 0.0904636982907597 | 7.99716717457639e-10 | 18.5439569215409 | 0 0 0 0
 chr1 | 11612 | 98 | 0.4489796 | 0.48979592 | 0.30612245 | 20 | 11.397959 | 10.5 | 12.333333 | 312 | 280 | 0.527027 | 1.3151919 | 0.0562953395595597 | 8.33420025699269e-16 | 31.8438504003201 | 0 0 0 0


### .clusters.txt
### .cov.gw.txt


## Example files


### Run forward simulations with SeDuS


### Process Output files
