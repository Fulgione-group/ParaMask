#!/bin/bash
java -jar /PATH_TO_PARAMASK_BYTECODE/ParaMask_Cluster_Seeds.jar\
	--cov PATH_TO_COV_FILE\
	--het PATH_TO_HET_FILE\
	--covgw PATH_TO_COVGW_FILE\
	--cutoff  DISTANCE_CUTOFF\
	--range CHR_START,CHR_END
