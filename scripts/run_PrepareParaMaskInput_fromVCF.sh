#!/bin/bash
java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/PrepareParaMaskInput_fromVCF.jar\
        --vcf $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf\
        --missingness 0.1
