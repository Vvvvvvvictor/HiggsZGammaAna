#!/bin/bash
#pogs=("BTV" "JME" "EGM" "MUO" "TAU")
pogs=("BTV" "EGM" "MUO" "TAU")
eras=("2016preVFP" "2016postVFP" "2017" "2018")
period="Run2"
campaign="UL"

for pog in "${pogs[@]}";
do

    if [[ ${pog} == "BTV" ]]; then jsonnames=("bjets"); fi
    if [[ ${pog} == "EGM" ]]; then jsonnames=("electron" "photon"); fi
    if [[ ${pog} == "MUO" ]]; then jsonnames=("muon"); fi
    if [[ ${pog} == "TAU" ]]; then jsonnames=("tau"); fi
    
    for era in "${eras[@]}";
    do

	for jsonname in "${jsonnames[@]}";
	do
	    correction --html ${pog}_${jsonname}_${era}_${campaign}.html summary ../../jsonpog-integration/POG/${pog}/${era}_${campaign}/${jsonname}.json
	    dir="/Volumes/dfs/Websites/c/cms-nanoAOD-integration/commonJSONSFs/"${pog}_${jsonname}_${period}_${campaign}
	    if [ ! -d "${dir}" ]; then mkdir ${dir}; fi
	    mv ${pog}_${jsonname}_${era}_${campaign}.html ${dir}
	    rm ${dir}/._*
	done
	
    done
    
done

