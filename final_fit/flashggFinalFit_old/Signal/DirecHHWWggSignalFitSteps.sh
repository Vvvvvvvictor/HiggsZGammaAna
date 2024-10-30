#!/bin/bash

## Abe Tishelman-Charny 
## 3 November 2019 
## Create HHWWgg Signal Model, datacard and run combine from HHWWgg Tagger output files located at $signalDirectory

# Need input file name format: X<mass>_HHWWgg_qqlnu.root

# Example usage:

# . DirecHHWWggSignalFitSteps.sh -i /eos/user/a/atishelm/ntuples/HHWWgg/testoutput/ -r HHWWgg_Signals -f
# This takes /eos/user/a/atishelm/ntuples/HHWWgg/testoutput/ as an input directory, and names all corresponding outputs with RedoX250Signal, and runs only the signalfit step (-k). 

cmsenv 

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o i:r:fskpdc -l help,inputDirectory:,procs:,smears:,massList:,scales:,scalesCorr:,useSSF:,useDCB_1G:,scalesGlobal:,flashggCats:,ext:,fTestOnly,calcPhoSystOnly,sigFitOnly,sigPlotsOnly,intLumi:,batch: -- "$@")
then
# something went wrong, getopt will put out an error message for us
echo "should exit"
fi
set -- $options

#signalFile=""
signalDirectory=""
runName=""
dofTest="false"
doPhoSyst="false"
createPoints="false"
doSignalFit="false"
plotSignal="false"
makeDatacard="false"
runCombine="false"

while [ $# -gt 0 ]
do
case $1 in
# -h|--help) usage; exit 0;;
#-i|--inputFile) signalFile=$2; shift ;;
-i|--inputDirectory) signalDirectory=$2; shift ;; 
-r|--run) runName=$2; shift ;;
-f|--fTest) dofTest="true"; shift ;;
-s|--phosys) doPhoSyst="true"; shift ;;
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

fggfinalfitDirec="/afs/cern.ch/work/a/atishelm/8Octflashggfinalfit/CMSSW_7_4_7/src/flashggFinalFit/"
sigDir=$fggfinalfitDirec
sigDir+="Signal"

# signalFile=$1
# runName=$2

#fileDir=${signalFile%/*}
fileDir=$signalDirectory
#fileDir+="/"
#echo "signalFile = $signalFile"
#echo "fileDir = $fileDir"

# steps: ftest, 120 and 130 gev points, Signalfit, plot, datacard, combine 
# do each step for each file in directory

##--FTest 
if [ $dofTest == 'true' ]
then
echo "signalDirectory: $signalDirectory"
signalDirectory+="*"
for f in $signalDirectory
do
  echo "Processing $f file..."
  signalFile=$f

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi

    #runName=${f::-5} # remove .root from what will be created folder name
    #runName="$(cut -d'/' -f-1 <<<$f)" 
    runName=${f##*/} # This trims everything from the front until a '/', greedily.
    runName=${runName::-5}
    fTestOutput=$fggfinalfitDirec
    fTestOutput+="Signal/"
    fTestOutput+=$runName
    datfilename=$fggfinalfitDirec
    datfilename+="Signal/dat/"
    datfilename+=$runName
    datfilename+=".dat"
    mass="$(cut -d'_' -f1 <<<$runName)"
    HHWWggLabel="${mass}_WWgg_qqlnugg"
    # HHWWggLabel=$runName
    # HHWWggLabel+="gg"
    # HHWWggLabel=X250_W #Wgg_qqlnugg
    echo "signalFile: $signalFile"
    echo "fTestOutput: $fTestOutput"
    echo "datfilename: $datfilename"
    echo "HHWWggLabel: $HHWWggLabel"
    ./bin/signalFTest -i $signalFile -p ggF --HHWWggLabel $HHWWggLabel -f HHWWggTag_0 -o $fTestOutput --datfilename $datfilename

done
fi

##-- Produce photon systematics dat file 
if [ $doPhoSyst == 'true' ]
then
  signalDirectory+="*"
  for f in $signalDirectory
  do
  echo "Processing $f file..."
  signalFile=$f

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi
  runName=${f##*/} # This trims everything from the front until a '/', greedily.
  runName=${runName::-5}
  mass="$(cut -d'_' -f1 <<<$runName)"
  HHWWggLabel="${mass}_WWgg_qqlnugg"
  outFile="${HHWWggLabel}.dat"

  echo "signalFile: $signalFile"
  echo "outFile: $outFile"
  echo "HHWWggLabel: $HHWWggLabel"

  ./bin/calcPhotonSystConsts -i $signalFile -o $outFile -p ggF -s HighR9EB,HighR9EE,LowR9EB,LowR9EE -r HighR9EBPhi,HighR9EBRho,HighR9EEPhi,HighR9EERho,LowR9EBPhi,LowR9EBRho,LowR9EEPhi,LowR9EERho -S MaterialCentralBarrel,MaterialForward,FNUFEB,ShowerShapeHighR9EE,ShowerShapeHighR9EB,ShowerShapeLowR9EE,ShowerShapeLowR9EB  -D outdir_HHWWgg_test -f HHWWggTag_0 -v 1 --HHWWggLabel $HHWWggLabel

  done 
fi
#./bin/calcPhotonSystConsts -i <input files comma separated> -o <output file> -p <comma separated process names> -s <photon energy scale categories> -r <photon energy scale categories>  -D <output dir for plots> -f <comma separated tag names>


# echo "doSignalFit = $doSignalFit"

##-- Create 120 and 130 GeV points 
if [ $doSignalFit == 'true' ]
then
signalDirectory+="*"
for f in $signalDirectory
do
  echo "Processing $f file..."

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi

    signalFile=$f
    runName=${f##*/} # This trims everything from the front until a '/', greedily.
    runName=${runName::-5}
    mass="$(cut -d'_' -f1 <<<$runName)" 
    HHWWggLabel="${mass}_WWgg_qqlnugg"
    
    python DirecShiftHiggsDatasets.py $fileDir $mass $HHWWggLabel

    sigFiles=""
    for m in 120 125 130
    do
        sigFiles+="${fileDir}X_signal_${mass}_${m}_HHWWgg_qqlnu.root,"
        # sigFiles+=$fileDir
        # sigFiles+="X_signal_"
        # sigFiles+=$mass
        # sigFiles+="_"
        # sigFiles+=$m
        # sigFiles+="_HHWWgg_qqlnu.root,"
    done        

    sigFiles=${sigFiles: : -1} # remove extra ","
    outputFile="CMS-HGG_sigfit_${mass}.root"
    # outputFile="CMS-HGG_sigfit"
    # outputFile+="_"
    # outputFile+=$mass
    # outputFile+=".root"

    phosysdat="${HHWWggLabel}.dat"

    echo "outputFile: $outputFile"
    echo "mass: $mass"
    echo "sigFiles: $sigFiles"
    echo "runName.dat: $runName.dat"
    echo "phosysdat: $phosysdat"

    ##-- SignalFit 
    ./bin/SignalFit -i $sigFiles -o $outputFile -p ggF -f HHWWggTag_0 --HHWWggLabel $HHWWggLabel -d dat/$runName.dat -s $phosysdat --procs ggF -v 1 --changeIntLumi 1 # key changeintlumi to 1 here 
  

#./bin/SignalFit -i $sigFiles -p ggF -f SL -d dat/$runName.dat -s empty.dat --procs ggF --changeIntLumi 1 # one inverse femtobarn     
# ./bin/SignalFit -i $sigFiles -p ggF -f SL -d dat/$runName.dat -s empty.dat --procs ggF --changeIntLumi 42.17
done
fi

##-- Plot 
if [ $plotSignal == 'true' ]
then
signalDirectory+="*"
for f in $signalDirectory
do

  echo "Processing $f file..."

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi
    runName=""
    runName=${f##*/} # This trims everything from the front until a '/', greedily.
    runName=${runName::-5}
    
mass="$(cut -d'_' -f1 <<<$runName)"

inputFile=""
inputFile+="CMS-HGG_sigfit_"
inputFile+=$mass
inputFile+=".root"

echo "runName: $runName"
echo "inputFile: $inputFile"
  #./bin/makeParametricSignalModelPlots -i CMS-HGG_sigfit.root  -o $runName -p ggF -f SL
	./bin/makeParametricSignalModelPlots -i $inputFile  -o $runName -p ggF -f HHWWggTag_0
done
fi 

###-- Datacard 
if [ $makeDatacard == 'true' ]
then

origDirec=$fggfinalfitDirec
origDirec+="Signal"

signalDirectory+="*"
for f in $signalDirectory
do

  echo "Processing $f file..."

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi

    runName=""
    runName=${f##*/} # This trims everything from the front until a '/', greedily.
    runName=${runName::-5}
    
    mass="$(cut -d'_' -f1 <<<$runName)"

    echo "mass: $mass"
    echo "runName: $runName"
    #echo "inputFile: $inputFile"

    datacardDir=$fggfinalfitDirec
    datacardDir+="Datacard"

    cd $datacardDir

    inputSignal=$fggfinalfitDirec
    inputSignal+="Signal/CMS-HGG_sigfit_"
    inputSignal+=$mass
    inputSignal+=".root"
    
    HHWWggLabel="${mass}_WWgg_qqlnugg"
    phosysdat="${HHWWggLabel}.dat"

    datacardName=$runName
    datacardName+=".dat"
    photonCatScales="${fggfinalfitDirec}Signal/${phosysdat}"

    echo "photonCatScales: $photonCatScales"

    python test_makeParametricModelDatacardFLASHgg.py -i $inputSignal -o $datacardName -p ggF -c HHWWggTag_0 --photonCatScales $photonCatScales --isMultiPdf --intLumi 41.5

    cd $origDirec 

# python test_makeParametricModelDatacardFLASHgg.py -i $inputSignal -o $datacardName -p ggF -c SL --photonCatScales $photonCatScales --isMultiPdf
done
fi

# ##-- Combine 
if [ $runCombine == 'true' ]
then
    # combineDir=$fggfinalfitDirec
    # combineDir+="Plots/FinalResults"
combineDir="/afs/cern.ch/work/a/atishelm/4NovCombineUpdated/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit"
origDirec=$fggfinalfitDirec
origDirec+="Signal"
signalDirectory+="*"
plotDirectory="$fggfinalfitDirec"
plotDirectory+="/Plots/FinalResults"
for f in $signalDirectory
do

  echo "Processing $f file..."

        if [[ $f == *"signal"* ]]; then
                echo "file contains 'signal', skipping"
                continue
        fi

    runName=""
    runName=${f##*/} # This trims everything from the front until a '/', greedily.
    runName=${runName::-5}

    mass="$(cut -d'_' -f1 <<<$runName)"

    echo "mass: $mass"

    datacardName=$runName
    datacardName+=".dat" 

    cd $combineDir
    cmsenv
    # cp ${fggfinalfitDirec}Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root

    cp ${fggfinalfitDirec}Signal/CMS-HGG_sigfit_${mass}.root CMS-HGG_sigfit_data_ggF_HHWWggTag_0.root
    cp ${fggfinalfitDirec}Background/HHWWgg_Background.root CMS-HGG_mva_13TeV_multipdf.root
    cp ${fggfinalfitDirec}Datacard/$datacardName CMS-HGG_mva_13TeV_datacard.txt

    # cp ../../Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root
    # cp ../../Signal/CMS-HGG_sigfit.root CMS-HGG_sigfit_data_ggF_SL.root
    # cp ../../Background/HHWWgg_Background.root CMS-HGG_mva_13TeV_multipdf.root 
    # cp ../../Datacard/$datacardName CMS-HGG_mva_13TeV_datacard.txt

    # combine  -M Asymptotic -m 125.00 --cminDefaultMinimizerType=Minuit2 -L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so  -d CMS-HGG_mva_13TeV_datacard.txt --run=blind -v 2
    combine CMS-HGG_mva_13TeV_datacard.txt -m 125 -M AsymptoticLimits --run=blind -v 2

    ### Copy combine output to plotting area 
    # configure plotter to take all .root files.
    # create script to run all of these fggrunjobs final steps should be easy like 10 lines  

    # combine  -M Asymptotic -m 125.00 --cminDefaultMinimizerType=Minuit2 -L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisGBRLikelihood.so  -d CMS-HGG_mva_13TeV_datacard.txt --run=expected  -v 2
    # What's the difference between --run=expected and --run=blind ?

    outName=$plotDirectory
    outName+="/higgsCombine.AsymptoticLimits.mH125."
    outName+=$mass
    outName+=".root"

    echo "combine output file: $outName"

    cp higgsCombineTest.AsymptoticLimits.mH125.root $outName
    
    cd $origDirec
    # cd $sigDir
    done
fi 
