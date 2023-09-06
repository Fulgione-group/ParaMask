#!/bin/bash
/netscratch/dep_coupland/grp_fulgione/bastiaan/software/SeDuS_1.10_exe Sim_SC_b$2_l$1_c0.05\
 -s 1\
 -z 100\
 -k 1000\
 -n 10000\
 -b $1\
 -u 0.001\
 -r 10\
 -c 0.05
