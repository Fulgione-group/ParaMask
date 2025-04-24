#!/bin/bash
Rscript --vanilla $PATH_TO_INSTALLATION_FOLDER/inst/scripts/run_ParaMask_EM.R\
      	--het $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.het.stat.txt\
        --missingness 0.3\
	--outdir $PATH_TO_INSTALLATION_FOLDER/Example_files/Output\
	--ID testing
