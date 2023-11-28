# Plot

## Codes for HZGamma plots
**1. Setup environment**
```
conda env create -f environment.yml -n hzgPlots
conda deactivate
conda activate hzgPlots

git clone -b HZGamma_run3 https://github.com/zebingwang/Plot.git
```

**2. Prepare Input Datasets**
Crate the input file```-i ./Data/2017/datasets.txt```, which has two rows (Datasets name and Datasets path). The python code ```scripts/PrepareDatasets.py``` will convert the parquet files into root files, which has the output path ```-o ./Data/2017```.

```
python scripts/PrepareDatasets.py -i ./Data/2017/datasets.txt -o ./Data/2017 --hadd
```

**3. Run the Plot Codes**
The main code is ```scripts/HZGamma_plot_Categorization.py```, which has all the needed libs in ```./lib```. ```./lib/Plot_Configs.py``` contains the definations of lumi, variables you want to plot, and the plot style.```./lib/Analyzer_Configs.py.py``` contains the input datasets information. ```./lib/Plot_Helper.py``` contains all the functions you need in the main plot codes.

```
python scripts/HZGamma_plot_Categorization.py -y [dataset year] --region [default:full] [--ln] [--ele]
```

## HZGamma Z Constrain Refit

**1. Prepare input Gen dataset**
Prepare the input ggH/VBFH Gen datasets, which includes the generator level informations about H->Zg->llg. You need to run the HiggsDNA and keep the generator level informations without any cuts. Once you got the ```.parquet``` file, you could write the path into ```.Data/2017/ZReFit/datasets.txt```. Then run:

```
python scripts/PrepareDatasets.py -i ./Data/2017/ZReFit/datasets.txt -o ./Data/2017/ZReFit/data_genZ
```

, which will convert ```.parquet``` to root files for the next processing.

**2. Get the Z true lineshape for the fit**

Get the fitted Z true lineshape using CB+3Gauss functions. You need the specify the input root file, which contain the generator level informations (produced by the first step), using parameter ```--genFile```. You should also specify the output path using parameter ```-o```, the channels using ```-c```. The default fitting is unbin fit, you could also use parameter ```--binningFit``` to speed up the fit.

```
python scripts/HZGamma_ZReFit.py --trueLineshapeFit --genFile ./Data/2017/ZReFit/data_genZ/ggH.root -o ./Data/2017/ZReFit/ -c ele/mu  (--binningFit --nBins 100)
```

**3. Do the Z Constrain Refit**
You could use ```--verbose``` to remove all the output information, but be sure there is no bug in your fit before you use it! Parameter ```--sampleName``` specified which sample you want to fit,. The default setting is ggH.

```
python scripts/HZGamma_ZReFit.py --doReFit --truelineshape ./Data/2017/ZReFit
```

You could also split the job by using ```--nJobs``` and ```--iJob```. After that, you could submit these jobs into Condor to speed up the refit. But, remember to modify the conda init to your own env in ```scripts/runbatch_ZreFit.jdl```.

```
condor_submit scripts/runbatch_ZreFit.jdl
```
Once all the jobs finished, hadd all the jobs.
```
hadd -f ggH.root ggH_job*.root
```


**4. Plot the Z Constrain Refit Results**

```
python scripts/HZGamma_ZReFit.py --plotReFit -c ele/mu
```