# HiggsZGammaAna
This is the full analysis chain of Run3 Higgs to ZGamma analysis. It is a combination of HiggsDNA, machine learning, spurial test and final fit. 

```
git clone https://github.com/Vvvvvvvictor/HiggsZGammaAna.git --recursive
```
## HiggsDNA for tagging
**1. Setup environment**
The conda env would occupy a large quota. Make sure that you have a few GB left in your path(exp. : `/eos/user`). This step would take a few hours, please keep patient.
```
cd HiggsZGammaAna/HiggsDNA/
conda env create -f environment.yml -p <some_path_where_you_have_more_disk_space>/.conda/envs/
```
Then, activate the conda environment

```
conda activate higgs-zg-ana
```

You may also want to increase your disk quota at [this link](https://resources.web.cern.ch/resources/Manage/EOS/Default.aspx), otherwise you may run out of space while installing your `conda` environment.

One additional package, `correctionlib`, must be installed via `pip`, rather than `conda`. Run
```
source setup.sh
```
to install this script.

**2. Add HiggsDNA package**
Install it by:
```
pip install -e .
```
If you notice issues with the `conda pack` command for creating the tarball, try updating and cleaning your environment with (after running `conda activate higgs-dna`):
```
conda env update --file environment.yml --prune
```

**3. Some useful note**
If you have some question of useage or code structure, please look at [HiggsDNA contents](https://sam-may.github.io/higgs_dna_tutorial.github.io/).

The name list of dataset can find in DAS system [example](https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FVBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer19UL17NanoAODv2-106X_mc2017_realistic_v8-v1%2FNANOAODSIM). Then you can put them in `metadata/samples/zgamma_tutorial.json`.

If you are comfused with the cross-section, please surf this [link](https://xsdb-temp.app.cern.ch/?columns=38289920&currentPage=0&ordDirection=1&ordFieldName=process_name&pageSize=50). If you want to find the path of dataset, please use this [link](https://cms-pdmv.cern.ch/mcm/requests?page=0&mcdb_id=102X_dataRun2_v11&shown=127).

## Machine learning for categorization
**1. Enter analysis environment**
Put this environment in HiggsDNA environment to reduce the space occupied by this program.
```
cd HiggsZGammaAna/hzgml
conda activate higgs-zg-ana
```
**2. Install packages**
Install the packages for machine learning(BDT and DNN).
```
pip install -r requirement.txt
```

**3. Usage of codes in /hzgml/**


## Spurial signel test
**1. Setup environment**
Get the CMS env first. Please make sure you have enter `HiggsZGammaAna' directory.
```
scram project CMSSW CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```
There are another two package needed. And need to clone them from github.
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit -b 112x-comb2021
cd HiggsAnalysis/CombinedLimit/
cd ../../
git clone https://github.com/cms-analysis/CombineHarvester.git HiggsAnalysis/CombineHarvester
cd HiggsAnalysis/CombineHarvester/
git checkout v2.0.0
cd ../../
scram b -j 9
```
Please note that we need to initialize **each time we setup a terminal**.
```
cd CMSSW_11_3_4/src
cmsenv
```
