ALP Flashgg Final Fit
========================
Contacts:

Repositories:
------------

Cloning the Repository
---------------
```
export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_13

cd CMSSW_10_2_13/src

cmsenv

git cms-init
```

Install the GBRLikelihood package which contains the RooDoubleCBFast implementation
```
git clone git@github.com:jonathon-langford/HiggsAnalysis.git
```
Install Combine as per the documentation here: cms-analysis.github.io/HiggsAnalysis-CombinedLimit/
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
```
Compile external libraries
-----------------------
```
cd ../HiggsAnalysis

cmsenv

scramv1 b clean; scramv1 b
```
Install Flashgg Final Fit packages
-----------------------
```
cd ..

git clone https://github.com/zebingwang/flashggFinalFit.git

cd flashggFinalFit/
```

Compile Signal Model Pakage
--------------

```
cd Signal

cmsenv

make clean

make

```

Compile Background Model Pakage
--------------

```
cd Background

cmsenv

make clean

make

```

Prepare dataset
--------------

Copy all the input datasets into dir: ```InputData```


Run
-----

```
./fit_HZGamma_model_UL.sh
```
