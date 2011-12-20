#!/bin/bash

cd /scratch/blanchou

#tar -xf CelegansAnalysis.tar
#tar -xvf TmpData.tar

#cd CelegansAnalysis

#./Cluster/mymex MEX/bilinear_mex.c
#./Cluster/mymex MEX/dp_score_mex.c
#./Cluster/mymex MEX/gaussian_mex.c MEX/gaussian_smooth.c
#./Cluster/mymex MEX/median_mex.c MEX/ctmf.c
#./Cluster/mymex MEX/imadm_mex.c MEX/gaussian_smooth.c

#cd ..

#mv CelegansAnalysis/Cluster/install_RECOS.m .
#mv OK-GZ920_MERGE_RECOS-.mat GZ920_MERGE_RECOS-.mat

#if [ ! -d TmpData ]
#  then mkdir TmpData
#fi

#mv tmpmat*.tmp TmpData

#if [ -d .matlab ]
#  then rm -r .matlab
#fi

#mkdir -m777 .matlab
#mkdir -m777 .matlab/R2009A
#mkdir -m777 .matlab/R2009b
#mkdir -m777 .matlab/R2009B
#mkdir -m777 .matlab/R2008A
#mkdir -m777 .matlab/R2008B

#pwd
#ls -lhR .matlab
#ls -lahR
eval "$MATLAB -nodisplay -nosplash -r ml_domain\($CPUUID\)"

