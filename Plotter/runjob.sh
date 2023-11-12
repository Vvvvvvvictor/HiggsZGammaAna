#!/bin/bash
/bin/hostname
gcc -v
pwd
export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh


#### job

##########hep_sub runjob.sh -g cms -mem 8000 -wt mid -o job.out -e job.err
#python ALP_plot_param.py -m -y 2016 

#python ALP_Optimization.py -y run2 -o ./optimize_run2 --doOpt -c 5
#python ALP_plot_param.py -y run2 -m --ln

#python ALP_plot_param.py -y run2 -m -S #--ln #--cut --mA M30
#python ALP_plot_param.py -y run2 -m 
#python ALP_plot_bkgCorr.py -y run2 -m
if [ $1 -eq 0 ]; then
    #python ALP_plot_param.py -y run2 -m --ele --ln  #--cut --cutVal 0.5 --mA M30 -b
    python ALP_plot_param.py -y run2 -m --ln
elif [ $1 -eq 1 ]; then
    #python ALP_plot_param.py -y run2 -m --mu --ln #--cut --cutVal 0.5 --mA M30 -b
    python ALP_plot_param.py -y run2  -m #--cut --cutVal 0.5 --mA M30 -b
elif [ $1 -eq 2 ]; then
    #python ALP_plot_param.py -y run2 -m --ele -C 
    python ALP_plot_param.py -y run2 -m -C --ln
elif [ $1 -eq 3 ]; then
    #python ALP_plot_param.py -y run2 -m --mu -C 
    python ALP_plot_param.py -y run2 -m -C
elif [ $1 -eq 4 ]; then
    #python ALP_plot_param.py -y run2 -m --ele -S 
    python ALP_plot_param.py -y run2 -m -S --ln
elif [ $1 -eq 5 ]; then
    #python ALP_plot_param.py -y run2 -m --mu -S 
    python ALP_plot_param.py -y run2 -m -S 
elif [ $1 -eq 6 ]; then
    #python ALP_Optimization.py -y run2 -o ./optimize_run2UL -p --sigVSscore -s --ele --doOpt -c 1
    python ALP_plot_bkgCorr.py -y run2 -m
elif [ $1 -eq 7 ]; then
    python ALP_Optimization.py -y run2 -o ./optimize_run2UL -p --sigVSscore -s --mu --doOpt -c 1
    python ALP_plot_bkgCorr.py -y run2 -m --mu
elif [ $1 -eq 8 ]; then
    #python ALP_Optimization.py -y run2 -o ./optimize_run2UL_ele -p --sigVSscore -s --ele --doOpt -c 2
    python ALP_BDTSys.py -y $2 -m --mu
elif [ $1 -eq 9 ]; then
    #python ALP_Optimization.py -y run2 -o ./optimize_run2UL_mu -p --sigVSscore -s --mu --doOpt -c 2
    python ALP_BDTSys.py -y $2 -m --ele
elif [ $1 -eq 10 ]; then
    python ALP_NormalizationSys_param.py -y run2 -m --ele
    python ALP_NormalizationSys_param.py -y $2 -m --mu
elif [ $1 -eq 11 ]; then
    python ALP_NormalizationSys_param.py -y $2 -m --ele
    #python ALP_BDTSys.py -y run2 --mu
elif [ $1 -eq 12 ]; then
    python ALP_plot_param.py -y run2  -m --cut --mA M10 --ln
elif [ $1 -eq 13 ]; then
    python ALP_plot_param.py -y run2  -m --cut --mA M20 --ln
elif [ $1 -eq 14 ]; then
    python ALP_plot_param.py -y run2  -m --cut --mA M30 --ln
fi