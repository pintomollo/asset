.\tar.exe -xvf CelegansAnalysis.tar
.\tar.exe -xvf TmpData.tar

cd CelegansAnalysis

.\mymex.bat MEX\bilinear_mex.c
.\mymex.bat MEX\dp_score_mex.c
.\mymex.bat MEX\gaussian_mex.c MEX\gaussian_smooth.c
.\mymex.bat MEX\median_mex.c MEX\ctmf.c
.\mymex.bat MEX\imadm_mex.c MEX\gaussian_smooth.c

cd ..

"%MATLAB%" -logfile output.txt -nojvm -nosplash -minimize -r install_RECOS "%CPUUID%"
