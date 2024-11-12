baseDir=`pwd`

cd $baseDir/Trees2WS/
#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # 设置输入目录
# inputDir="/eos/user/j/jiehan/root/output_cor_syst"
# outputWSDir="/eos/user/j/jiehan/root/ws_cor_syst"

# # 遍历文件夹
# for proc in "$inputDir"/*; do
#     for year in "$proc"/*; do
#         for file in "$year"/*.root; do
#             if [[ -f "$file" ]]; then
#                 # 获取年份和文件名
#                 proc_name=$(basename "$proc")
#                 year_num=$(basename "$year")
#                 # 生成命令
#                 cmd="python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile $file --productionMode $proc_name --year $year_num --outputWSDir $outputWSDir --doSystematics" 
#                 echo $cmd
#                 $cmd
#             fi
#         done
#     done
# done

# python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile "/eos/user/j/jiehan/root/output_cor_syst/ggH/2017/output.root" --productionMode ggH --year 2017 --outputWSDir $outputWSDir --doSystematics

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# inputDir="/eos/user/j/jiehan/root/output_bkg"
# outputWSDir="/eos/user/j/jiehan/root/ws_data"

# for proc in "$inputDir"/Data*; do
#     for year in "$proc"/*; do
#         for file in "$year"/*.root; do
#             if [[ -f "$file" ]]; then
#                 # 获取年份和文件名
#                 proc_name=$(basename "$proc")
#                 year_num=$(basename "$year")
#                 # 生成命令
#                 cmd="python3 trees2ws_data.py --inputConfig config_bkg.py --inputTreeFile $file --outputWSDir $outputWSDir" 
#                 echo $cmd
#                 $cmd
#             fi
#         done
#     done
# done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Signal/
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode fTest
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode calcPhotonSyst
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode getDiagProc
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat # using n gaussian functions
# rm -rf outdir_2017/signalFit_nGauss; mv outdir_2017/signalFit outdir_2017/signalFit_nGauss
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--useDCB"
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
#     python3 RunPlotter.py --years 2017 --ext 2017 --cats $cat --groupSignalFitByProc
# done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Background/

cmsenv
# # ~~~~One Time Setup~~~~
make clean
make -j 16
# # ~~~~~~~~~End~~~~~~~~~~

# # sleep 1
# # clear

rm -rf outdir_2017
python3 RunBackgroundScripts.py --inputConfig config_2017.py --mode fTestParallel

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Datacard/
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
#     python3 makeYields.py --inputWSDirMap 2017=/eos/user/j/jiehan/root/ws_cor_syst/signal_2017/  --inputBkgWSDirMap 2017=/eos/user/j/jiehan/root/ws_data/Data_2017/ --sigModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_2017/signalFit/output --bkgModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/outdir_2017 --year 2017 --cat $cat
#     python3 makeDatacard.py --years 2017 --cat $cat
# done

# python3 makeYields.py --inputWSDirMap 2017=/eos/user/j/jiehan/root/ws_cor_syst/signal_2017/  --inputBkgWSDirMap 2017=/eos/user/j/jiehan/root/ws_data/Data_2017/ --sigModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_2017/signalFit/output --bkgModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/outdir_2017 --year 2017 --cat ggH0 -procs 'ggH,VBF'
# python3 makeDatacard.py --years 2017 -cat ggH0

# Combine the datacards

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Combine/
# rm -rf output_Datacard
# cp -r $baseDir/Datacard/output_Datacard/ ./
# rm -rf output_plots
# if [ ! -d "output_plots" ]; then
#     mkdir output_plots
# fi

# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
#     python3 RunText2Workspace.py --year 2017 --cat $cat
# done

# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
#     cd output_plots
#     combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d ../output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --setParameterRanges r=-4,6 --points 101 --firstPoint 0 --lastPoint 100 -n _scan1D_syst_r_${cat}
#     cd ../
#     plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_${cat}.MultiDimFit.mH125.root --y-cut 8 --y-max 8 --output output_plots/r_statsyst_${cat} --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" #--others output_plots/higgsCombine_scan1D_syst_r_${cat}.MultiDimFit.mH125.root:"Stat only":2
# done


# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 lep VH ZH ttHh ttHl;do
#     combine output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root -M Significance --toysFrequentist -t -1 --expectSignal=1 -m 125.0 -n $cat 
#     # combine output_Datacard/Datacard_2017_mu_inclusive.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat
#     # python3 RunFits.py --year 2017 --inputWS output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --inputJson inputs.json --mode mu_inclusive
# done

# python3 RunText2Workspace.py -year 2017
# combine output_datacard_rootfile_/Datacard_2017_mu_inclusive.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n allCats --toysFrequentist
# combine output_datacard_rootfile_/Datacard_2017_mu_inclusive.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n allCats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Plots/
# python makeSplusBModelPlot.py --inputWSFile ../Combine/output_datacard_rootfile_/Datacard_2017_mu_inclusive.root --cats all 
