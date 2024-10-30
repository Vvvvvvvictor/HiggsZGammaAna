#!/bin/bash

year=( 16 17 18 )
nYear=${#year[@]}
massList=( 5 15 30 )
nMass=${#massList[@]}

cd ./Signal/ALP_SigModel_param/
mkdir fit_results_runII
cd fit_results_runII

for ((iBin=0; iBin<$nMass; iBin++))
    do
    mkdir M${massList[$iBin]}
    cd M${massList[$iBin]}

    for ((jBin=0; jBin<$nYear; jBin++))
        do
        #pwd
        #echo "../../fit_results_${year[$jBin]}/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt ./datacard_ALPmass${massList[$iBin]}_20${year[$jBin]}.txt"
        #echo "../../fit_results_${year[$jBin]}/M${massList[$iBin]}/CMS-HGG_mva_13TeV_multipdf.root ./CMS-HGG_mva_13TeV_multipdf_20${year[$jBin]}.root"
        #echo "../../fit_results_${year[$jBin]}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0.root ./CMS-HGG_sigfit_data_ggh_cat0_20${year[$jBin]}.root"
        #echo '#####'

        cp ../../fit_results_${year[$jBin]}/M${massList[$iBin]}/datacard_ALPmass${massList[$iBin]}.txt ./datacard_ALPmass${massList[$iBin]}_20${year[$jBin]}.txt
        cp ../../fit_results_${year[$jBin]}/M${massList[$iBin]}/CMS-HGG_mva_13TeV_multipdf.root ./CMS-HGG_mva_13TeV_multipdf_20${year[$jBin]}.root
        cp ../../fit_results_${year[$jBin]}/M${massList[$iBin]}/CMS-HGG_sigfit_data_ggh_cat0.root ./CMS-HGG_sigfit_data_ggh_cat0_20${year[$jBin]}.root
        sed -i "s/CMS-HGG_mva_13TeV_multipdf.root/CMS-HGG_mva_13TeV_multipdf_20${year[$jBin]}.root/g" ./datacard_ALPmass${massList[$iBin]}_20${year[$jBin]}.txt
        sed -i "s/CMS-HGG_sigfit_data_ggh_cat0.root/CMS-HGG_sigfit_data_ggh_cat0_20${year[$jBin]}.root/g" ./datacard_ALPmass${massList[$iBin]}_20${year[$jBin]}.txt

    done
    combineCards.py datacard_ALPmass${massList[$iBin]}_2016.txt datacard_ALPmass${massList[$iBin]}_2017.txt datacard_ALPmass${massList[$iBin]}_2018.txt &> datacard_ALPmass${massList[$iBin]}.txt
    combine datacard_ALPmass${massList[$iBin]}.txt -M AsymptoticLimits --run=blind -m 125 --rAbsAcc 0.00000001
    cd ../
done

mkdir plot
cd plot
cp ../../../ALP_plot/*.py ./
sed -i "15s/ALP_SigModel/ALP_SigModel_param/" com_plot.py
sed -i "15s/fit_results_16/fit_results_runII/" com_plot.py
sed -i "34s/.*/lumi_sqrtS = '35.9 + 41.5 + 56.9 = 134.3 fb^{-1}'/" CMS_lumi.py
python com_plot.py -fb
