#!/bin/bash
export OPENBLAS_NUM_THREADS=40 OMP_NUM_THREADS=40 MKL_NUM_THREADS=40
Rscript --vanilla ~/PATH_TO_PARAMASK_SRC/ParaMask_EM_v2.4.R "--hetfile" PATH_TO_HET_FILE \
"--outpath" PATH_TO_OUTPUTDIR "--missingness" MISSINGNESS "--verbose" "--ID" RUN_ID
