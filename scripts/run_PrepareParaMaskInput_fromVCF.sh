#!/bin/bash
java -jar /PATH_TO_PARAMASK_BYTECODE/PrepareParaMaskInput_fromVCF.jar\
	-VF\
	--vcf INPUT_VCF.sh\
	--missingness 0.1

