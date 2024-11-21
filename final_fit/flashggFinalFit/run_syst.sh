baseDir=`pwd`

cd $baseDir/Trees2WS/
#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 设置输入目录
inputDir="/eos/user/j/jiehan/root/output_cor_syst"
outputWSDir="/eos/user/j/jiehan/root/ws_cor_syst"

# 遍历文件夹
for proc in "$inputDir"/*; do
    for year in "$proc"/*; do
        for file in "$year"/*.root; do
            if [[ -f "$file" ]]; then
                # 获取年份和文件名
                proc_name=$(basename "$proc")
                year_num=$(basename "$year")
                # 生成命令
                cmd="python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile $file --productionMode $proc_name --year $year_num --outputWSDir $outputWSDir --doSystematics" 
                echo $cmd
                $cmd
            fi
        done
    done
done

# python3 trees2ws.py --inputConfig config_signal.py --inputTreeFile "/eos/user/j/jiehan/root/output_cor_syst/ggH/2017/output.root" --productionMode ggH --year 2017 --outputWSDir $outputWSDir --doSystematics

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inputDir="/eos/user/j/jiehan/root/output_bkg"
outputWSDir="/eos/user/j/jiehan/root/ws_data"

for proc in "$inputDir"/Data*; do
    for year in "$proc"/*; do
        for file in "$year"/*.root; do
            if [[ -f "$file" ]]; then
                # 获取年份和文件名
                proc_name=$(basename "$proc")
                year_num=$(basename "$year")
                # 生成命令
                cmd="python3 trees2ws_data.py --inputConfig config_bkg.py --inputTreeFile $file --outputWSDir $outputWSDir" 
                echo $cmd
                $cmd
            fi
        done
    done
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Signal/
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode fTest
python3 RunSignalScripts.py --inputConfig config_2017.py --mode calcPhotonSyst
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode getDiagProc
# python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat # using n gaussian functions
# rm -rf outdir_2017/signalFit_nGauss; mv outdir_2017/signalFit outdir_2017/signalFit_nGauss
python3 RunSignalScripts.py --inputConfig config_2017.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--useDCB"
cp -r outdir_2017/signalFit outdir_2017/signalFit_DCB
for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
    python3 RunPlotter.py --years 2017 --ext 2017 --cats $cat --groupSignalFitByProc
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Background/

cmsenv
# # ~~~~One Time Setup~~~~
make clean
make -j 16
# # ~~~~~~~~~End~~~~~~~~~~

# # sleep 1
# # clear

# rm -rf outdir_2017/bkgfTest-Data/
python3 RunBackgroundScripts.py --inputConfig config_2017.py --mode fTestParallel

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Datacard/

for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
    python3 makeYields.py --inputWSDirMap 2017=/eos/user/j/jiehan/root/ws_cor_syst/signal_2017/  --inputBkgWSDirMap 2017=/eos/user/j/jiehan/root/ws_data/Data_2017/ --sigModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_2017/signalFit/output --bkgModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/outdir_2017 --year 2017 --cat $cat --doSystematics
    python3 makeDatacard.py --years 2017 --cat $cat --doSystematics
done

# python3 makeYields.py --inputWSDirMap 2017=/eos/user/j/jiehan/root/ws_cor_syst/signal_2017/  --inputBkgWSDirMap 2017=/eos/user/j/jiehan/root/ws_data/Data_2017/ --sigModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/outdir_2017/signalFit/output --bkgModelWSDir /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/outdir_2017 --year 2017 --cat ggH0 -procs 'ggH,VBF'
# python3 makeDatacard.py --years 2017 -cat ggH0

# Combine the datacards
combineCards.py ggH0=output_Datacard/pruned_datacard_2017_ggH0_.txt ggH1=output_Datacard/pruned_datacard_2017_ggH1_.txt ggH2=output_Datacard/pruned_datacard_2017_ggH2_.txt ggH3=output_Datacard/pruned_datacard_2017_ggH3_.txt > output_Datacard/pruned_datacard_2017_ggH_.txt
combineCards.py VBF0=output_Datacard/pruned_datacard_2017_VBF0_.txt VBF1=output_Datacard/pruned_datacard_2017_VBF1_.txt VBF2=output_Datacard/pruned_datacard_2017_VBF2_.txt VBF3=output_Datacard/pruned_datacard_2017_VBF3_.txt > output_Datacard/pruned_datacard_2017_VBF_.txt # VBF2=output_Datacard/pruned_datacard_2017_VBF2_.txt VBF3=output_Datacard/pruned_datacard_2017_VBF3_.txt
combineCards.py ggH0=output_Datacard/pruned_datacard_2017_ggH0_.txt ggH1=output_Datacard/pruned_datacard_2017_ggH1_.txt ggH2=output_Datacard/pruned_datacard_2017_ggH2_.txt ggH3=output_Datacard/pruned_datacard_2017_ggH3_.txt VBF0=output_Datacard/pruned_datacard_2017_VBF0_.txt VBF1=output_Datacard/pruned_datacard_2017_VBF1_.txt VBF2=output_Datacard/pruned_datacard_2017_VBF2_.txt VBF3=output_Datacard/pruned_datacard_2017_VBF3_.txt VH=output_Datacard/pruned_datacard_2017_VH_.txt ZH=output_Datacard/pruned_datacard_2017_ZH_.txt ttHh=output_Datacard/pruned_datacard_2017_ttHh_.txt ttHl=output_Datacard/pruned_datacard_2017_ttHl_.txt > output_Datacard/pruned_datacard_2017_.txt # ggH0=output_Datacard/pruned_datacard_2017_ggH0_.txt ggH1=output_Datacard/pruned_datacard_2017_ggH1_.txt ggH2=output_Datacard/pruned_datacard_2017_ggH2_.txt VH=output_Datacard/pruned_datacard_2017_VH_.txt ZH=output_Datacard/pruned_datacard_2017_ZH_.txt ttHh=output_Datacard/pruned_datacard_2017_ttHh_.txt ttHl=output_Datacard/pruned_datacard_2017_ttHl_.txt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Combine/
rm -rf output_Datacard
cp -r $baseDir/Datacard/output_Datacard/ ./
rm -rf output_plots
if [ ! -d "output_plots" ]; then
    mkdir output_plots
fi

# ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl

for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
    python3 RunText2Workspace.py --year 2017 --cat $cat --common_opts "-m 125 higgsMassRange=122,128 --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints cms"
done

# ----------------------------------mu scan----------------------------------

for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
    cd output_plots

    # combine -M MultiDimFit --mass 125 ../output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --saveWorkspace -n ${cat}

    combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d ../output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_${cat} --freezeParameters MH
    combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d ../output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_${cat}.statOnly --freezeParameters allConstrainedNuisances,MH

    cd ../
    plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_${cat}.MultiDimFit.mH125.root --y-cut 8 --y-max 8 --output output_plots/r_statsyst_${cat} --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" --others output_plots/higgsCombine_scan1D_syst_r_${cat}.statOnly.MultiDimFit.mH125.root:"Stat. Only":2 #output_plots/higgsCombine_scan1D_syst_r_${cat}_p0to6.MultiDimFit.mH125.root:"r0to6":4
done

# combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d comb_hzg.root --setParameterRanges r=-4,6 --points 101 -n allcat

# combine --algo grid -M MultiDimFit -m 125.0 -d higgsCombineVBF0.MultiDimFit.mH125.root --setParameterRanges r=-4,6 --points 10 -n _scan1D_syst_r_VBF0

# plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_VBF0.MultiDimFit.mH125.root --y-cut 8 --y-max 8 --output output_plots/r_statsyst_VBF0 --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" --others output_plots/higgsCombine_scan1D_syst_r_VBF0_p-4to0.MultiDimFit.mH125.root:"r-4to0":2 output_plots/higgsCombine_scan1D_syst_r_VBF0_p0to6.MultiDimFit.mH125.root:"r0to6":4

# ---------------------------For ggH, VBF and all categories---------------------------
# # text2workspace.py -m 125.38 --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints cms inputs/12538/comb_2021_hzg_lhcomb.txt -o comb_hzg.root
text2workspace.py output_Datacard/pruned_datacard_2017_ggH_.txt -o output_datacard_rootfile_/Datacard_2017_ggH_.root -m 125 higgsMassRange=122,128 --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints cms
text2workspace.py output_Datacard/pruned_datacard_2017_VBF_.txt -o output_datacard_rootfile_/Datacard_2017_VBF_.root -m 125 higgsMassRange=122,128 --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints cms
text2workspace.py output_Datacard/pruned_datacard_2017_.txt -o output_datacard_rootfile_/Datacard_2017_.root -m 125 higgsMassRange=122,128 --for-fits --no-wrappers --X-pack-asympows --optimize-simpdf-constraints cms

cd output_plots
combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d ../output_datacard_rootfile_/Datacard_2017_ggH_.root --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_ggH.statOnly --freezeParameters allConstrainedNuisances
combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d ../output_datacard_rootfile_/Datacard_2017_ggH_.root --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_ggH --freezeParameters MH #-v 3 > log_ggH.log

# combine -M MultiDimFit --mass 125 ../output_datacard_rootfile_/Datacard_2017_VBF_.root --saveWorkspace -n VBF --floatOtherPOIs 0 --expectSignal 1 -t -1 --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2
combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_VBF_.root --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_VBF.statOnly --freezeParameters allConstrainedNuisances 
combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 --setParameterRanges r=0,10 --points 21 -n _scan1D_syst_r_VBF --freezeParameters MH  -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_VBF_.root #-v 3 > log_nonuis_01.log
#-d $baseDir/Combine/output_plots/higgsCombineVBF.MultiDimFit.mH125.root --snapshotName MultiDimFit --skipInitialFit

combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100 --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_.root --setParameterRanges r=0,3 --points 51 -n _scan1D_syst_r_statOnly --freezeParameters allConstrainedNuisances,MH #,pdfindex_ggH0_13TeV,pdfindex_ggH1_13TeV,pdfindex_ggH2_13TeV,pdfindex_ggH3_13TeV,pdfindex_VBF0_13TeV,pdfindex_VBF1_13TeV,pdfindex_VBF2_13TeV,pdfindex_VBF3_13TeV,pdfindex_VH_13TeV,pdfindex_ZH_13TeV,pdfindex_ttHh_13TeV,pdfindex_ttHl_13TeV
combine --floatOtherPOIs 0 --expectSignal 1 -t -1 -P r --algo grid --alignEdges 1 --saveSpecifiedNuis all --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100 --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 -M MultiDimFit -m 125.0 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_.root --setParameterRanges r=0,3 --points 51 -n _scan1D_syst_r_ --freezeParameters MH #,pdfindex_ggH0_13TeV,pdfindex_ggH1_13TeV,pdfindex_ggH2_13TeV,pdfindex_ggH3_13TeV,pdfindex_VBF0_13TeV,pdfindex_VBF1_13TeV,pdfindex_VBF2_13TeV,pdfindex_VBF3_13TeV,pdfindex_VH_13TeV,pdfindex_ZH_13TeV,pdfindex_ttHh_13TeV,pdfindex_ttHl_13TeV
# pdfindex_ggH0_13TeV,pdfindex_ggH1_13TeV,pdfindex_ggH2_13TeV,pdfindex_ggH3_13TeV,pdfindex_VBF0_13TeV,pdfindex_VBF1_13TeV,pdfindex_VBF2_13TeV,pdfindex_VBF3_13TeV,pdfindex_VH_13TeV,pdfindex_ZH_13TeV,pdfindex_ttHh_13TeV,pdfindex_ttHl_13TeV
# ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl
cd ../

# plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_ggH.MultiDimFit.mH125.root --y-cut 8 --y-max 8 --output output_plots/r_statsyst_ggH --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" --others output_plots/higgsCombine_scan1D_syst_r_ggH.statOnly.MultiDimFit.mH125.root:"Stat. Only":2
# plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_VBF.MultiDimFit.mH125.root --y-cut 8 --y-max 8 --output output_plots/r_statsyst_VBF --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" --others output_plots/higgsCombine_scan1D_syst_r_VBF.statOnly.MultiDimFit.mH125.root:"Stat. Only":2
# plot1DScan.py output_plots/higgsCombine_scan1D_syst_r_.MultiDimFit.mH125.root --y-cut 4 --y-max 4 --output output_plots/r_statsyst_allCats --POI r --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --logo-sub "Preliminary" --others output_plots/higgsCombine_scan1D_syst_r_statOnly.MultiDimFit.mH125.root:"Stat. Only":2
plot1DScan.py /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_plots/more_env_func/higgsCombine_scan1D_syst_r_statOnly.MultiDimFit.mH125.root --y-cut 4 --y-max 4 --output output_plots/r_statsyst_allCats_compare --POI r --translate ../Plots/pois_mu.json --main-label "More Func.(fix index)" --main-color 1 --logo-sub "Preliminary" --others /afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/final_fit/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_plots/one_env_func/higgsCombine_scan1D_syst_r_statOnly.MultiDimFit.mH125.root:"One Func. in Envelope":2

# ---------------------------------For Significance---------------------------------

for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
    combine output_Datacard/pruned_datacard_2017_${cat}_.txt -M Significance --toysFrequentist -t -1 --expectSignal=1 -m 125.0 -n $cat 
    # combine output_Datacard/Datacard_2017_mu_inclusive.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n $cat
    # python3 RunFits.py --year 2017 --inputWS output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --inputJson inputs.json --mode mu_inclusive
done

combine output_Datacard/pruned_datacard_2017_ggH_.txt -M Significance --toysFrequentist -t -1 --expectSignal=1 -m 125.0 -n ggH
combine output_Datacard/pruned_datacard_2017_VBF_.txt -M Significance --toysFrequentist -t -1 --expectSignal=1 -m 125.0 -n VBF
combine output_Datacard/pruned_datacard_2017_.txt -M Significance --toysFrequentist -t -1 --expectSignal=1 -m 125.0 -n allCats

# # python3 RunText2Workspace.py -year 2017
# # combine output_datacard_rootfile_/Datacard_2017_mu_inclusive.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n allCats --toysFrequentist
# # combine output_datacard_rootfile_/Datacard_2017_mu_inclusive.root -M AsymptoticLimits --run=blind -m 125.0 --rAbsAcc 0.00000001 -n allCats

# ---------------------------Impact of systematics on mu---------------------------
# if [ ! -d "impacts" ]; then
#     mkdir impacts
# fi
# cd impacts
# # for cat in ggH3 VBF0 VBF1 VBF2 VBF3 VH ZH ttHh ttHl;do
# #     combineTool.py -M Impacts -m 125 --doInitialFit --robustFit 1 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --parallel 10 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2
# #     combineTool.py -M Impacts -m 125 --robustFit 1 --doFits -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root --parallel 10
# #     combineTool.py -M Impacts -m 125 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_${cat}__mu_inclusive.root -o impacts_${cat}.json
# #     plotImpacts.py -i impacts_${cat}.json -o impacts_${cat}
# # done
# # cd ..

# combineTool.py -M Impacts -m 125 --doInitialFit --robustFit 1 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_.root --parallel 10
# combineTool.py -M Impacts -m 125 --robustFit 1 --doFits -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_.root --parallel 10
# combineTool.py -M Impacts -m 125 -d $baseDir/Combine/output_datacard_rootfile_/Datacard_2017_.root -o impacts_allCats.json
# plotImpacts.py -i impacts_allCats.json -o impacts_allCats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split Line~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd $baseDir/Plots/
python makeSplusBModelPlot.py --inputWSFile ../Combine/output_datacard_rootfile_/Datacard_2017_.root --cats all 
