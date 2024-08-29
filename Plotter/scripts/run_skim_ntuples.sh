#!/bin/bash                                                                                                                                                                       
echo "==============STARTED=============="

#i_low=$(( 100*($1 - 1) + 1 ))
#i_hig=$(( 100*($1) ))

#for ((i=$i_low; i<=$i_hig; i++))
#    do
#    python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/data_run2_2016postVFP_CR/Data_2016postVFP/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/data_CR_16postVFP/data_CR_16postVFP_${i}.root


#done

if [ $1 -eq 1 ];then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_2016postVFP_CR/DYJetsToLL_2016postVFP/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2016postVFP.root
    for ((i=1; i<=50; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_17/DY_2017_${i}.root
    done
elif [ $1 -eq 2 ]; then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_2016preVFP_CR/DYJetsToLL_2016preVFP/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2016preVFP.root
    for ((i=51; i<=100; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_17/DY_2017_${i}.root
    done
elif [ $1 -eq 3 ];then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_2016postVFP_CR/DYJetsToLL_2016postVFP/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2016postVFP.root
    for ((i=101; i<=155; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_17/DY_2017_${i}.root
    done
elif [ $1 -eq 4 ]; then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2017.root
    for ((i=1; i<=50; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_18_CR/DYJetsToLL_2018/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_18/DY_2018_${i}.root
    done
elif [ $1 -eq 5 ]; then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2017.root
    for ((i=51; i<=100; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_18_CR/DYJetsToLL_2018/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_18/DY_2018_${i}.root
    done
elif [ $1 -eq 6 ]; then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_17_CR/DYJetsToLL_2017/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2017.root
    for ((i=101; i<=150; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_18_CR/DYJetsToLL_2018/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_18/DY_2018_${i}.root
    done
elif [ $1 -eq 7 ]; then
    #python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_18_CR/DYJetsToLL_2018/merged_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_data_Normalizing/DY_2018.root
    for ((i=151; i<=205; i++))
        do
        python /afs/cern.ch/work/z/zewang/private/HZGamma/HiggsZGammaAna/hzgml/scripts/skim_ntuples.py -i /afs/cern.ch/work/z/zewang/private/HZGamma/bkgmc_DY_18_CR/DYJetsToLL_2018/job_${i}/output_job_${i}_nominal.parquet -o /eos/user/z/zewang/HZGamma_data/run2UL_CR/DY_18/DY_2018_${i}.root
    done
fi



echo "==============FINISHED==========="
