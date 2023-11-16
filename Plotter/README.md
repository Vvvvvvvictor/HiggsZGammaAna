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
