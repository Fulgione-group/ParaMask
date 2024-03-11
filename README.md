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
  Chromosome	Position	Non.missing	Minor.allele.freq	Heterozygous.geno.freq	Homozygous1.geno.freq	Homozygous2.geno.freq	Mean.coverage	Mean.coverage.hom	Mean.coverage.het	Het.reference.allele.depth	Het.alt.allele.depth	Het.allele.ratio  Het.allele.deviation	L1	L2	LLR	EM_class	allele.deviation.seed	finalClass	cluster"
  chr1	10641	100	0.025	0.05	0.95	0	19.68	19.68421	19.6	72	26	0.7346939	4.646702	0.411559811962468	0.963941228136795	-0.851075965471138	2	1  1	1
  chr1	10680	100	0.04	0.08	0.92	0	19.63	19.554348	20.5	118	46	0.7195122	5.6222553	0.22685182553148	0.943842011392948	-1.42566173794247	2	1  1	1
  chr1	10776	99	0.08080808	0.14141414	0.8484849	1	20.818182	20.964706	19.928572	185	94	0.6630824	5.4480276	0.203133317411145	0.0113927604667986	2.88088439259977	2	1	1	1
  chr1	10798	100	0.49	0.94	0.02	4	20.9	14	21.340425	993	1013	0.49501497	-0.4465443	2.90395537401198e-14	0.0210172995916636	-27.3077081818396	2	0  1	1
  chr1	10811	99	0.05050505	0.1010101	0.8989899	0	19.959597	19.853933	20.9	143	66	0.68421054	5.3262014	0.147898455696237	0.931023048035841	-1.83975810506303	2	1	1	1
  chr1	10906	100	0.015	0.03	0.97	0	19.87	19.927835	18	38	16	0.7037037	2.993821	0.596228304224322	0.977961587648625	-0.494846738335215	2	1  1	1
  chr1	10971	99	0.5	1	0	0	21.424242	0	21.424242	1063	1058	0.5011787	0.10856746	2.3785162706989e-19	0.608498643720493	-42.3858792938336	20	1	1
  chr1	11065	100	0.025	0.05	0.95	0	19.96	19.83158	22.4	80	32	0.71428573	4.5355735	0.411559811962468	0.963941228136795	-0.851075965471138	2	1  1	1
  chr1	11226	99	0.26767677	0.35353535	0.5555556	9	10.686869	10.453125	11.114285	193	196	0.49614397	-0.15210603	0.0904636982907597	7.99716717457639e-10	18.5439569215409	0	0	0	0
  chr1	11612	98	0.4489796	0.48979592	0.30612245	20	11.397959	10.5	12.333333	312	280	0.527027	1.3151919	0.0562953395595597	8.33420025699269e-16	31.8438504003201	0	0	0	0
  chr1	11993	99	0.43434343	0.4040404	0.36363637	23	9.474748	9.338983	9.675	185	202	0.47803617	-0.8641586	0.0423081470021819	6.12813943557946e-19	38.773449969909	0	0	0	0
  chr1	12373	99	0.030303031	0.060606062	0.93939394	0	10.565657	10.634409	9.5	30	27	0.5263158	0.3973597	0.338789194072882	0.957120910950128	-1.03855165964169	1	0	0	0
  chr1	12397	98	0.010204081	0.020408163	0.97959185	0	10.897959	10.895833	11	9	13	0.4090909	-0.8528029	0.711835410020391	0.985169275283179	-0.324966760499041	1	0	0	0
  chr1	12415	95	0.12105263	0.17894737	0.7894737	3	9.263158	9.051282	10.235294	74	100	0.42528737	-1.9710549	0.157212930371547	1.73152084290818e-05	9.11377119459786	0	0	0	0
  chr1	12578	99	0.005050505	0.01010101	0.989899	0	10.30303	10.295918	11	3	8	0.27272728	-1.5075567	0.845957500000001	0.992514	-0.159771996317597	1	0	0	0
  chr1	12677	95	0.005263158	0.010526316	0.9894737	0	10.3052635	10.308511	10	5	5	0.5	0	0.845827900000008	0.992514	-0.15992520723932	1  0	0  0
  chr1	12738	98	0.0051020407	0.010204081	0.9897959	0	10.193877	10.092784	20	10	10	0.5	0	0.845926100000007	0.992514	-0.15980911471037210	0	0
  chr1	12789	100	0.46	0.36	0.36	28	10.15	9.859375	10.666667	196	188	0.5104167	0.4082483	0.0164352534466298	9.90104677390362e-23	46.5584900017412	0  0	0	0
  chr1	12871	98	0.44897962	0.5510204	0.1734694	27	9.836735	9.386364	10.203704	267	284	0.4845735	-0.7242243	0.027760645469523	3.00467028127724e-13	25.2493024784942	0	0	0	0
  chr1	13600	96	0.13020833	0.17708333	0.78125	4	10.0625	9.797468	11.294118	97	95	0.5052083	0.14433756	0.104148040117188	7.83616779655431e-07	11.7974038082202	0	0	0	0
  chr1	13628	96	0.015625	0.03125	0.96875	0	10.40625	10.397849	10.666667	21	11	0.65625	1.767767	0.595413477781447	0.977961587648625	-0.496214307989407	1	0	0	0
  chr1	13965	97	0.005154639	0.010309278	0.9896907	0	10.835052	10.8125	13	6	7	0.46153846	-0.2773501	0.845894100000007	0.992514	-0.159846943789377	1	0	0	0
  ```



### .clusters.txt
### .cov.gw.txt


## Example files


### Run forward simulations with SeDuS


### Process Output files
