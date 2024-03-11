#!/bin/bash
Rscript --vanilla $PATH_TO_INSTALLATION_FOLDER/ParaMask/ParaMask_EM_v0.2.5.R\
      	--het $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.het.stat.txt\
        --missingness 0.1\
	--outdir $PATH_TO_INSTALLATION_FOLDER/Example_files/Output\
	--ID test
