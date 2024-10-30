#!/bin/bash

## Abe Tishelman-Charny 
## 3 November 2019 
## Create HHWWgg Signal Model, datacard and combine from flashgg EventDumper output file located at $signalFile


# Need input file name format: X<mass>_HHWWgg_qqlnu.root

# Example usage:

# . HHWWgg_Signal_Fit_Steps.sh -i /eos/user/a/atishelm/ntuples/HHWWgg/HHWWgg_v1_Dumper_Signal_AllEvents_WithWorkspace/X250_HHWWgg_qqlnu.root -r RedoX250Signal -k
# This takes /eos/.../X250_HHWWgg_qqlnu.root as an input file, and names all corresponding outputs with RedoX250Signal, and runs only the signalfit step (-k). 

cmsenv 

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o i:r:fkpdc -l help,inputFile:,procs:,smears:,massList:,scales:,scalesCorr:,useSSF:,useDCB_1G:,scalesGlobal:,flashggCats:,ext:,fTestOnly,calcPhoSystOnly,sigFitOnly,sigPlotsOnly,intLumi:,batch: -- "$@")
then
# something went wrong, getopt will put out an error message for us
echo "should exit"
fi
set -- $options

signalFile=""
runName=""
dofTest="false"
createPoints="false"
doSignalFit="false"
plotSignal="false"
makeDatacard="false"
runCombine="false"

while [ $# -gt 0 ]
do
case $1 in
# -h|--help) usage; exit 0;;
-i|--inputFile) signalFile=$2; shift ;;
-r|--run) runName=$2; shift ;;
-f|--fTest) dofTest="true"; shift ;;
-k|--signalfit) doSignalFit="true"; shift ;; 
-p|--plotsig) plotSignal="true"; shift ;;
-d|--dcard) makeDatacard="true"; shift ;;
-c|--runComb) runCombine="true"; shift ;;

(--) shift; break;;
(-*) usage; echo "$0: [ERROR] - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
(*) break;;
esac
shift
done

fggfinalfitDirec="/afs/cern.ch/work/c/chuw/HHWWgg/flashggFinalFit/CMSSW_10_2_13/src/flashggFinalFit/"
#fggfinalfitDirec="/afs/cern.ch/work/a/atishelm/8Octflashggfinalfit/CMSSW_7_4_7/src/flashggFinalFit/"
sigDir=$fggfinalfitDirec
sigDir+="Signal"

# signalFile=$1
# runName=$2
datfilename=$fggfinalfitDirec
datfilename+="Signal/dat/"
datfilename+=$runName
datfilename+=".dat"

fileDir=${signalFile%/*}
fileDir+="/"
echo "signalFile = $signalFile"
echo "fileDir = $fileDir"

# steps: ftest, 120 and 130 gev points, Signalfit, plot, datacard, combine 

##--FTest 
if [ $dofTest == 'true' ]
then
    fTestOutput=$fggfinalfitDirec
    fTestOutput+="Signal/"
    fTestOutput+=$runName
    ./bin/signalFTest -i $signalFile -p ggF -f SL -o $fTestOutput --datfilename $datfilename
fi

##-- Produce photon systematics dat file 
# Doing stat only for now 

# echo "doSignalFit = $doSignalFit"

##-- Create 120 and 130 GeV points 
if [ $doSignalFit == 'true' ]
then
    #python shiftHiggsDatasets.py $fileDir 

    sigFiles=""
    for m in 120 125 130
    do
        sigFiles+=$fileDir
        sigFiles+="X_signal_250_"
        sigFiles+=$m
        sigFiles+="_HHWWgg_qqlnu.root,"
    done        

    sigFiles=${sigFiles: : -1} # remove extra ","
    echo "$sigFiles"
    ##-- SignalFit 
    # ./bin/SignalFit -i $sigFiles -p ggF -f SL -d dat/$runName.dat -s empty.dat --procs ggF 
    ./bin/SignalFit -i $sigFiles -p ggF -f HHWWggTag_0 -d dat/$runName.dat -s empty.dat --procs ggF --changeIntLumi 1 # one femtobarn 
    # ./bin/SignalFit -i $sigFiles -p ggF -f SL -d dat/$runName.dat -s empty.dat --procs ggF --changeIntLumi 42.17
fi

##-- Plot 
if [ $plotSignal == 'true' ]
then
    ./bin/makeParametricSignalModelPlots -i CMS-HGG_sigfit.root  -o $runName -p ggF -f SL
fi 

# ##-- Datacard 
if [ $makeDatacard == 'true' ]
then
    datacardDir=$fggfinalfitDirec
    datacardDir+="Datacard"

    cd $datacardDir

    inputSignal=$fggfinalfitDirec
    inputSignal+="Signal/CMS-HGG_sigfit.root"
    datacardName=$runName
    datacardName+=".dat"
    photonCatScales=$fggfinalfitDirec
    photonCatScales+="Signal/empty.dat"
    python test_makeParametricModelDatacardFLASHgg.py -i $inputSignal -o $datacardName -p ggF -c SL --photonCatScales $photonCatScales --isMultiPdf --intLumi 41.5
# python test_makeParametricModelDatacardFLASHgg.py -i $inputSignal -o $datacardName -p ggF -c SL --photonCatScales $photonCatScales --isMultiPdf
fi

# ##-- Combine 
if [ $runCombine == 'true' ]
then
    # combineDir=$fggfinalfitDirec
    # combineDir+="Plots/FinalResults"

    datacardName=$runName
    datacardName+=".dat" 

    combineDir="/afs/cern.ch/work/a/atishelm/4NovCombineUpdated/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit"

    cd $combineDir
    cmsenv
    cp ${fggfinalfitDirec}Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root
    cp ${fggfinalfitDirec}Background/HHWWgg_Background.root CMS-HGG_mva_13TeV_multipdf.root
    cp ${fggfinalfitDirec}Datacard/$datacardName CMS-HGG_mva_13TeV_datacard.txt
    # cp ../../Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root
    # cp ../../Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root
    # cp ../../Background/HHWWgg_Background.root CMS-HGG_mva_13TeV_multipdf.root 
    # cp ../../Datacard/$datacardName CMS-HGG_mva_13TeV_datacard.txt

    # combine  -M Asymptotic -m 125.00 --cminDefaultMinimizerType=Minuit2 -L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so  -d CMS-HGG_mva_13TeV_datacard.txt --run=blind -v 2
    combine CMS-HGG_mva_13TeV_datacard.txt -m 125 -M AsymptoticLimits --run=blind -v 2 

    # combine  -M Asymptotic -m 125.00 --cminDefaultMinimizerType=Minuit2 -L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so  -d CMS-HGG_mva_13TeV_datacard.txt --run=expected  -v 2
    # What's the difference between --run=expected and --run=blind ?

    # cd $sigDir
fi 
