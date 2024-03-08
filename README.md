# Overview
ParaMask consists of 3 steps
1. PrepareVCF_fromVCF.sh:  Computed 3 ParaMaskInput files from a VCF file
    - Het file: file with all relevant Information for the Expectation Maximization algorithm and Read Ratio Deviations
    - Cov file: Coverage file per sample and per site
    - CovGW file: Genome-wide coverage of non-missing sites per individual 
2. ParaMask_EM.sh: Utelizes population genomic signatures of multicopies for seed SNP generation.
    - EM algorithm to jointly fit 2 Beta-binomial regressions on the heterozygote frequency with varying Minor Allele frequency: One for single-copy regions, and one for multicopy regions
    - RRD testing: 
## Simulations

### Run forward simulations with SeDuS


### Process Output files
