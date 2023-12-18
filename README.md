---

[![Repobeats analytics image](https://repobeats.axiom.co/api/embed/bac34668655cc1118c86a0b1831cf095e159606f.svg "Repobeats analytics image")](https://github.com/Vvvvvvvictor/HiggsZGammaAna/pulse)

---
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
cd HiggsZGammaAna/HiggsDNA
conda deactivate
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

The logic computation sequance is: 1. `&`, 2. `|`. It would be better use `()` between each `&`.

If you have some question of useage or code structure, please look at [HiggsDNA contents](https://sam-may.github.io/higgs_dna_tutorial.github.io/).

The name list of dataset can find in DAS system [(example)](https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FVBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer19UL17NanoAODv2-106X_mc2017_realistic_v8-v1%2FNANOAODSIM). Then you can put them in `metadata/samples/zgamma_tutorial.json`.

If you are comfused with the cross-section of a MC dataset, please get it through [this link](https://xsdb-temp.app.cern.ch/?columns=38289920&currentPage=0&ordDirection=1&ordFieldName=process_name&pageSize=50). 

If you want to find the path of a dataset, please get it through [this link](https://cms-pdmv.cern.ch/mcm/requests?page=0&mcdb_id=102X_dataRun2_v11&shown=127).

If you want to get the golden json, these files, `/eos/user/c/cmsdqm/www/CAF/certification/Collisions*/Cert_Collisions*_*_*_Golden.json`, are recommended. Please put them in `metadata/golden_json/` and also change the golden json choice in the tagger you use.

You can find the name of variables in NanoAOD in [this link](https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv11/2022postEE/doc_Muon_Run2022G-PromptNanoAODv11_v1-v2.html)

## Z Constrain Refit

See details in ```Plotter```.

## Machine learning for categorization

### Setup environment
**1. Enter analysis environment**

Need another envirenment, you can set up it through these codes:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc9-opt/setup.sh
python3 -m venv /eos/(your dir)/hzgmlenv
```

You MUST make a new teminal window for the following steps, please change the `source /eos/(your dir)/hzgmlenv/bin/activate` in `scripts/install.sh`
```
cd HiggsZGammaAna/hzgml/
source scripts/install.sh
```

It would be better to put this environment in `/eos` or other disk but not `/afs`. You can set the path by changing it in the `install.sh` and `setup.sh` mentioned below. The important thing is that they are supposed to be the same in these two files.

**2. Set up environment each time**

You need to set up the environment each time you create a new terminal.
```
source scripts/setup.sh
```

### Scripts to run the tasks

#### Prepare the training inputs (skimmed ntuples)

The core script is `skim_ntuples.py`. The script will apply the given skimming cuts to input samples and produce corresponding skimmed ntuples. You can run it locally for any specified input files. In H->mumu, it is more convenient and time-efficient to use `submit_skim_ntuples.py` which will find all of the input files specified in `data/inputs_config.json` and run the core script by submitting the condor jobs.

- The script is hard-coded currently, meaning one needs to directly modify the core script to change the variable calculation and the skimming cuts.
- The output files (skimmed ntuples) will be saved to the folder named `skimmed_ntuples` by default.
- `submit_skim_ntuples.py` will only submit the jobs for those files that don't exist in the output folder.
- `run_skim_ntuples.sh` will run all skim jobs locally

```
sh scripts/run_skim_ntuples.sh
```
or
```
python scripts/skim_ntuples.py [-i input_file_path] [-o output_file_path]
```

#### Start XGBoost analysis!

The whole ML task consists of training, applying the weights, optimizing the BDT boundaries for categorization, and calculating the number counting significances. The wrapper script `run_all.sh` will run everything. Please have a look!

**1. Make some directories**
```
mkdir -p models outputs plots
```

**2. Training a model**

The training script `train_bdt.py` will train the model in four-fold, and transform the output scores such that the unweighted signal distribution is flat. The detailed settings, including the preselections, training variables, hyperparameters, etc, are specified in the config file `data/training_config.json`.

```
python scripts/train_bdt.py [-r TRAINED_MODEL] [-f FOLD] [--save]

Usage:
  -r, --region        The model to be trained. Choices: 'zero_jet', 'one_jet', 'two_jet' or 'VBF'.
  -f, --fold          Which fold to run. Default is -1 (run all folds)
  --save              To save the model into HDF5 files, and the pickle files
```

**3. Applying the weights**

Applying the trained model (as well as the score transformation) to the skimmed ntuples to get BDT scores for each event can be done by doing:
```
python scripts/apply_bdt.py [-r REGION]
```
The script will take the settings specified in the training config file `data/training_config.json` and the applying config file `data/apply_config.json`.

**4. Optimizing the BDT boundaries**

`categorization_1D.py` will take the Higgs classifier scores of the samples and optimize the boundaries that give the best combined significance. `categorization_2D.py`, on the other hand, takes both the Higgs classifier scores and the VBF classifier scores of the samples and optimizes the 2D boundaries that give the best combined significance.

```
python scripts/categorization_1D.py [-r REGION] [-f NUMBER OF FOLDS] [-b NUMBER OF CATEGORIES] [-n NSCAN] [--floatB] [--minN minN] [--skip]

Usage:
  -f, --fold          Number of folds of the categorization optimization. Default is 1.
  -b, --nbin          Number of BDT categories
  -n, --nscan         Number of scans. Default is 100
  --minN,             minN is the minimum number of events required in the mass window. The default is 5.
  --floatB            To float the last BDT boundary, which means to veto the lowest BDT score events
  --skip              To skip the hadd step (if you have already merged signal and background samples)
```

```
python scripts/categorization_2D.py [-r REGION] [-f NUMBER OF FOLDS] [-b NUMBER OF CATEGORIES] [-b NUMBER OF ggF CATEGORIES] [-n NSCAN] [--floatB] [--minN minN] [--skip]

Usage:
  -f, --fold          Number of folds of the categorization optimization. Default is 1.
  -b, --nbin          Number of BDT categories
  -n, --nscan         Number of scans. Default is 100
  --minN,             minN is the minimum number of events required in the mass window. The default is 5.
  --floatB            To float the last BDT boundary, which means to veto the lowest BDT score events
  --skip              To skip the hadd step (if you have already merged signal and background samples)
```

## Spurial signel test

For this step, there is no need to activate the conda analysis environment.

### Setup environment

Get the CMS env first. Please make sure you have enter `HiggsZGammaAna' directory.
```
cd HiggsZGammaAna/SSTest/
scram project CMSSW CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```
There are another two package needed. And need to clone them from github.
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit -b 112x-comb2021
git clone https://github.com/cms-analysis/CombineHarvester.git HiggsAnalysis/CombineHarvester -b v2.0.0
scram b -j 9
```
Please note that you need to initialize **each time you setup a terminal** by doing this.
```
cd CMSSW_11_3_4/src
cmsenv
```

### Scripts to run a task
**1. Prepare the background and signal template**

We can get a file(`bkg_sig_template.root`) with the bkg and sig shape in each channels and categories by running this.
```
python SSTest/Generate_template.py
```
You need to modify the range of the hist and make sure it is the same as that in `SSTest.cpp` or `SSTest_core_function.cpp`.

**2. Run spurial signal test and find the best bkg function**

Run this code:
```
root -l -q SSTest/SSTest.cpp
```
There is some option in `SSTest.cpp`, cat: which category of this channel to test, channel: one_jet and something similar, sig: How many times of signal injected, bkg_fun: function used to discribe background.

Also a Spurious signal test based on core function method is provided:
```
root -l -q SSTest/SSTest_core_function.cpp
```
There is some option in `SSTest_core_function.cpp`, cat: which category of this channel to test, channel: one_jet and something similar, sig: How many times of signal injected.

for the background function for fitting, we can change inside the codes.

**Take care that the observed variable in this code should be `CMS_hzg_mass`, because this is the name used in `ZGMCShape.root`. We should use the same name, or we would fail.**

We select the best background function by spurial signal, chi square and F test. We can see the fitting result and name of the best function in the log file.(not done yet)

## NOTE
After setting up three environment(HiggsDNA, machine learning and spurial signal test), we can set up those three by running the compacted setup shell script:
```
source setup.sh
```

## Final Fit
### Setup environment
**1.Setup CMSSW environment**
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
```

**2.Install dependency package**
```
git clone https://github.com/jonathon-langford/HiggsAnalysis.git
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
```

**3.Compile dependency package**
```
cd ../
scramv1 b clean; scramv1 b
```

**4.Install Flashgg Final Fit package**
```
cd final_fit/
mv flashggFinalFit CMSSW_10_2_13/src
```
setup CMSSW environment
```
cd CMSSW_10_2_13/src
cmsenv
cd flashggFinalFit/
```

```
cd Signal
make clean
make
```

```
cd Background
make clean
make
```
