#!/bin/bash
java -jar $PATH_TO_INSTALLATION_FOLDER/ParaMask/ParaMask_Cluster_Seeds.jar\
	--cov $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.cov.stat.txt\
	--het $PATH_TO_INSTALLATION_FOLDER/Example_files/Output/test_EMresults.het\
	--covgw  $PATH_TO_INSTALLATION_FOLDER/Example_files/Input/Simulations_10PercentDuplications.vcf.cov.gw.txt\
	--cutoff  $(tail -1 $PATH_TO_INSTALLATION_FOLDER/Example_files/Output/test_EMresults.dist)\
	--range 1,1000000
