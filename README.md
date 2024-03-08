# Overview
ParaMask consists of 3 steps
1. PrepareVCF_fromVCF.sh:  Computes 3 ParaMaskInput files from a VCF file
    - Het file: file with all relevant Information for the Expectation Maximization algorithm and Read Ratio Deviations
    - Cov file: Coverage file per sample and per site
    - CovGW file: Genome-wide coverage of non-missing sites per individual 
2. ParaMask_EM.sh: Utelizes population genomic signatures of multicopies for seed SNP generation.
    - EM algorithm to jointly fit 2 Beta-binomial regressions on the heterozygote frequency with varying Minor Allele frequency: One for single-copy regions, and one for multicopy regions.
       * Classifiaction based on the LLR
    - RRD testing: Takes the mean and variance from the read ratio deviation of single copy classified SNPs and constructs a normal confidence interval to test previously uncertain classified SNPs
    - EM algortihm to dissect distances of SNPs within and between multicopy regions.
       * Takes distances between seed SNPs and calculates mean parameters of a mixture of geometric_distributions.
       * We set the cutoff distance to where the the two geometrics have the same density.
       * To make this procedure more stable by default we repeat this procedure 1000 and take the median cutoff.
3. ParaMask_Cluster_Seeds.sh: Cluster SNPs into multicopy haplotypes and gives final annation of SNPs, using Seeds, coverage and distance. Computes three Outputfiles:
    - The original Het file with final status
    - Cluster file with each multicopy SNP and annotation
    - Bed file with single- and multicopy regions
## Simulations

### Run forward simulations with SeDuS


### Process Output files
