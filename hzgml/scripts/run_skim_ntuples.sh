#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

input="/eos/home-j/jiehan/parquet/nanov9/"
target="/eos/home-j/jiehan/root/skimmed_ntuples/"

for ((i=1; i<=27; i++))
    do
    python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_ZGToLLG_2016postVFP_CR/ZGToLLG_2016postVFP/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/ZG16postVFP/ZGToLLG_CR_16postVFP_${i}.root

done



echo "==============FINISHED==========="
