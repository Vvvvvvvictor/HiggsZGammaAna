#!/usr/bin/env python
import copy
import os
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd
import uproot
# from root_pandas import *
import pickle
from sklearn.preprocessing import StandardScaler, QuantileTransformer
import xgboost as xgb
from tqdm import tqdm
import logging
import gc  # Add garbage collection
from pdb import set_trace
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError + 1
pd.options.mode.chained_assignment = None

try:
    import psutil  # For memory monitoring
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    print("WARNING: psutil not available, memory monitoring disabled")

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    PYARROW_AVAILABLE = True
except ImportError:
    PYARROW_AVAILABLE = False
    print("WARNING: pyarrow not available, using pickle for temporary files")

from weighted_quantile_transformer import WeightedQuantileTransformer # New import

def get_memory_usage():
    """Get current memory usage in GB."""
    if PSUTIL_AVAILABLE:
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024 / 1024
    else:
        return 0.0

def check_hadd_available():
    """Check if ROOT's hadd command is available."""
    import subprocess
    try:
        result = subprocess.run(['which', 'hadd'], capture_output=True, text=True)
        return result.returncode == 0
    except:
        return False

HADD_AVAILABLE = check_hadd_available()

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', action='store', nargs=2, default=['data/training_config_BDT.json', 'data/apply_config_BDT.json'], help='Region to process')
    parser.add_argument('-i', '--inputFolder', action='store', default='/eos/home-j/jiehan/root/skimmed_ntuples_run2', help='directory of training inputs')
    parser.add_argument('-m', '--modelFolder', action='store', default='models', help='directory of BDT models')
    parser.add_argument('-o', '--outputFolder', action='store', default='/eos/home-j/jiehan/root/outputs/test', help='directory for outputs')
    parser.add_argument('-r', '--region', action='store', choices=['two_jet', 'one_jet', 'zero_jet', 'zero_to_one_jet', 'VH_ttH', 'all_jet'], default='zero_jet', help='Region to process')
    parser.add_argument('-cat', '--category', action='store', nargs='+', help='apply only for specific categories')
    parser.add_argument('-s', '--shield', action='store', type=int, default=-1, help='Which variables needs to be shielded')
    parser.add_argument('-a', '--add', action='store', type=int, default=-1, help='Which variables needs to be added')

    return parser.parse_args()

class ApplyXGBHandler(object):
    "Class for applying XGBoost"

    def __init__(self, configPath, region=''):

        print('===============================')
        print('  ApplyXGBHandler initialized')
        print('===============================')

        args=getArgs()
        self._shield = args.shield
        self._add = args.add

        self._region = region
        self._inputFolder = ''
        self._inputTree = region if region else 'inclusive'
        self._modelFolder = ''
        self._outputFolder = ''
        self._chunksize = 1000000
        self._category = []
        self._branches = []
        self._outbranches = []

        self.m_models = {}
        self.m_tsfs = {}

        self.train_variables = {}
        self.randomIndex = 'eventNumber'

        self.models = {}
        self.observables = []
        self.preselections = []

        self.readApplyConfig(configPath[1])
        self.readTrainConfig(configPath[0])
        self.arrangeBranches()
        # self.arrangePreselections()

    def readApplyConfig(self, configPath):
        """Read configuration file formated in json to extract information to fill TemplateMaker variables."""
        try:
            member_variables = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("_") and not attr.startswith('m_')]

            stream = open(configPath, 'r')
            configs = json.loads(stream.read())

            # read from the common settings
            config = configs["common"]
            for member in config.keys():
                if member in member_variables:
                    setattr(self, member, config[member])

            # read from the region specific settings
            if self._region:
                config = configs[self._region]
                for member in config.keys():
                    if member in member_variables:
                        setattr(self, member, config[member])
                if '+preselections' in config.keys():
                    self.preselections += config['+preselections']
                if '+observables' in config.keys():
                    self.observables += config['+observables']

        except Exception as e:
            logging.error("Error reading apply configuration '{config}'".format(config=configPath))
            logging.error(e)

    def readTrainConfig(self, configPath):

        try:
            stream = open(configPath, 'r')
            configs = json.loads(stream.read())
            # if (self._add>=0):
            #     configs["common"]["train_variables"].append(configs["common"]["+train_variables"][self._add])
   
            config = configs["common"]

            if 'randomIndex' in config.keys(): self.randomIndex = config['randomIndex']
 
            if self.models:
                for model in self.models:
    
                    # read from the common settings
                    config = configs["common"]
                    if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
    
                    # read from the region specific settings
                    if model in configs.keys():
                        config = configs[model]
                        if self._add >= 0:
                            config["+train_variables"].append(config["test_variables"][self._add])
                        if 'train_variables' in config.keys(): self.train_variables[model] = config['train_variables'][:]
                        if '+train_variables' in config.keys(): self.train_variables[model] += config['+train_variables']
                        if self._shield >= 0:
                            self.train_variables[model].pop(self._shield)

                        print("\n\n")
                        print(self.train_variables[model])
                        print(len(self.train_variables[model]))
                        print("\n\n")

        except Exception as e:
            logging.error("Error reading training configuration '{config}'".format(config=configPath))
            logging.error(e)

    def arrangeBranches(self):

        self._branches = set()
        for model in self.models:
            self._branches = self._branches | set(self.train_variables[model])

        self._branches = self._branches | set([self.randomIndex]) | set([p.split()[0] for p in self.preselections]) | set(self.observables)
        self._branches = list(self._branches)

        for model in self.models:
            self.train_variables[model] = [x.replace('noexpand:', '') for x in self.train_variables[model]]
        self.preselections = [x.replace('noexpand:', '') for x in self.preselections]
        self.randomIndex = self.randomIndex.replace('noexpand:', '')

        self._outbranches = [branch for branch in self._branches if 'noexpand' not in branch]

    def arrangePreselections(self):

        if self.preselections:
            self.preselections = ['data.' + p for p in self.preselections]

    def setInputFolder(self, inputFolder):
        self._inputFolder = inputFolder

    def setModelFolder(self, modelFolder):
        self._modelFolder = modelFolder

    def setOutputFolder(self, outputFolder):
        self._outputFolder = outputFolder

    def preselect(self, data):
        print(self.preselections)
        for p in self.preselections:
            data = data.query(p)

        return data

    def loadModels(self):

        if self.models:
            for model in self.models:
                print('XGB INFO: Loading BDT model: ', model)
                self.m_models[model] = []
                for i in range(4):
                    bst = xgb.Booster()
                    bst.load_model('%s/BDT_%s_%d.h5'%(self._modelFolder, model, i))
                    self.m_models[model].append(bst)
                    del bst

    def loadTransformer(self):
        
        if self.models:
            for model in self.models:
                print('XGB INFO: Loading score transformer for model: ', model)
                self.m_tsfs[model] = []
                for i in range(4):
                    tsf = pickle.load(open('%s/BDT_tsf_%s_%d.pkl'%(self._modelFolder, model, i), "rb" ), encoding = 'latin1' )
                    self.m_tsfs[model].append(tsf)

    def applyBDT(self, category, scale=1, shift=0):
        outputbraches = copy.deepcopy(self._outbranches)
        branches = copy.deepcopy(self._branches)
        # branches += ["Z_sublead_lepton_pt", "gamma_mvaID_WP80", "gamma_mvaID_WPL", "gamma_pt", "Z_pt"]
        outputbraches += []
        # if category == "DYJetsToLL":
        #     branches.append('n_iso_photons')
        # if category != "data_fake" and category != "mc_true" and category != "mc_med":
        #     branches.append('gamma_mvaID_WP80')
        # if category == "data_fake" or category == "mc_true" or category == "mc_med":
        #     branches += ['weight_err']
        #     outputbraches += ['weight_err']
        if category == "mc_true" or category == "mc_med":
            branches += ['tagger']
            outputbraches += ['tagger']

        # print(branches)
        # print(outputbraches)
        
        outputContainer = self._outputFolder + '/' + self._region
        output_path = outputContainer + '/%s.root' % category
        if not os.path.isdir(outputContainer): os.makedirs(outputContainer)
        if os.path.isfile(output_path): os.remove(output_path)

        f_list = []
        cat_folder = self._inputFolder + '/' + category
        for f in os.listdir(cat_folder):
            if f.endswith('.root'): f_list.append(cat_folder + '/' + f)

        print('-------------------------------------------------')
        for f in f_list: print('XGB INFO: Including sample: ', f)

        # Use ultra-low memory mode by default
        print("XGB INFO: Using ultra-low memory mode")
        self._process_files(f_list, output_path, branches, outputbraches, scale, shift, category)

    def _process_files(self, f_list, output_path, branches, outputbraches, scale, shift, category):
        """Process files with ultra-low memory usage."""
        import tempfile
        temp_files = []
        
        initial_memory = get_memory_usage()
        print(f"XGB INFO: Initial memory usage: {initial_memory:.2f} GB")
        
        # Process each file individually and save to temp files immediately
        for file_idx, filename in enumerate(tqdm(sorted(f_list), desc='XGB INFO: Applying BDTs to %s samples' % category, bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}')):
            try:
                file = uproot.open(filename)
            except Exception as e:
                print('XGB ERROR: Failed to open file: ', filename)
                continue
            
            # Process each chunk immediately without batching
            chunk_count = 0
            for data in file[self._inputTree].iterate(library='pd', step_size=self._chunksize):  # Use smaller chunks
                data = self.preselect(data)
                
                for i in range(4):
                    data_s = data[(data[self.randomIndex]-shift)%314159%4 == i]
                    if data_s.shape[0] == 0: continue
                    
                    data_o = data_s.copy()

                    for model in self.train_variables.keys():
                        x_Events = data_s[self.train_variables[model]]
                        dEvents = xgb.DMatrix(x_Events)
                        scores = self.m_models[model][i].predict(dEvents)
                        if len(scores) > 0:
                            scores_t = self.m_tsfs[model][i].transform(scores.reshape(-1,1)).reshape(-1)
                        else:
                            scores_t = scores
                    
                        xgb_basename = self.models[model]
                        data_o[xgb_basename] = scores
                        data_o[xgb_basename+'_t'] = scores_t
                    
                    # Save immediately to temp ROOT file
                    if len(data_o) > 0:
                        temp_file = tempfile.NamedTemporaryFile(suffix='.root', delete=False)
                        temp_file.close()
                        with uproot.recreate(temp_file.name) as root_file:
                            root_file[self._region] = data_o
                        temp_files.append(temp_file.name)
                    
                    # Clear immediately
                    del data_o
                    chunk_count += 1
                    
                    # Force garbage collection more frequently
                    if chunk_count % 4 == 0:
                        gc.collect()
                
                # Clear chunk data
                del data, data_s
                gc.collect()
            
            file.close()
            gc.collect()
            
            current_memory = get_memory_usage()
            print(f"XGB INFO: Processed file {file_idx+1}/{len(f_list)}, Memory: {current_memory:.2f} GB, Temp files: {len(temp_files)}")
        
        # Combine temp files with minimal memory usage
        print(f"XGB INFO: Combining {len(temp_files)} temporary files...")
        self._combine_temp_files(temp_files, output_path)
        
        # Clean up
        for temp_file in temp_files:
            try:
                os.remove(temp_file)
            except:
                pass
        
        final_memory = get_memory_usage()
        print(f"XGB INFO: Final memory usage: {final_memory:.2f} GB")
        
    def _combine_temp_files(self, temp_files, output_path):
        """Combine temporary files using ROOT hadd for optimal performance."""
        if not temp_files:
            return
        
        print(f"XGB INFO: Ultra-low memory combination of {len(temp_files)} files...")
        
        # Temporary files are already in ROOT format, no conversion needed
        print(f"XGB INFO: All {len(temp_files)} temporary files are already in ROOT format")
        
        # Use ROOT's hadd to combine all files
        if temp_files:
            if HADD_AVAILABLE:
                print("XGB INFO: Using ROOT hadd to combine files...")
                hadd_command = f"hadd -f {output_path} " + " ".join(temp_files)
                
                import subprocess
                try:
                    result = subprocess.run(hadd_command, shell=True, capture_output=True, text=True)
                    if result.returncode == 0:
                        print(f"XGB INFO: Successfully combined {len(temp_files)} files using hadd")
                    else:
                        print(f"XGB ERROR: hadd failed: {result.stderr}")
                        # Fallback to manual combination
                        self._fallback_combine_files(temp_files, output_path)
                except Exception as e:
                    print(f"XGB ERROR: Failed to run hadd: {e}")
                    # Fallback to manual combination
                    self._fallback_combine_files(temp_files, output_path)
            else:
                print("XGB WARNING: hadd not available, using manual combination...")
                self._fallback_combine_files(temp_files, output_path)
    
    def _fallback_combine_files(self, root_files, output_path):
        """Fallback method to combine ROOT files manually if hadd fails."""
        print("XGB INFO: Using fallback method to combine files...")
        
        if not root_files:
            return
        
        # Read first file to initialize
        with uproot.open(root_files[0]) as first_file:
            combined_data = first_file[self._region].arrays(library='pd')
        
        # Append remaining files one by one
        for root_file in root_files[1:]:
            try:
                with uproot.open(root_file) as file:
                    data = file[self._region].arrays(library='pd')
                    combined_data = pd.concat([combined_data, data], ignore_index=True, sort=False)
                    del data
                    gc.collect()
            except Exception as e:
                print(f"XGB WARNING: Failed to read {root_file}: {e}")
                continue
        
        # Write final result
        with uproot.recreate(output_path) as output_file:
            output_file[self._region] = combined_data
        
        del combined_data
        gc.collect()

    def _combine_files_streaming(self, temp_files, output_path):
        """Combine temporary files using ROOT hadd in streaming mode."""
        if not temp_files:
            return
        
        print(f"XGB INFO: Combining {len(temp_files)} files using ROOT hadd...")
        
        # Convert temporary files to ROOT files in batches
        root_temp_files = []
        batch_size = 20  # Process 20 temp files at a time
        
        for batch_start in range(0, len(temp_files), batch_size):
            batch_end = min(batch_start + batch_size, len(temp_files))
            batch_files = temp_files[batch_start:batch_end]
            
            # Create ROOT files for this batch
            batch_root_files = []
            for i, temp_file in enumerate(batch_files):
                try:
                    if PYARROW_AVAILABLE and temp_file.endswith('.parquet'):
                        df = pd.read_parquet(temp_file, engine='pyarrow')
                    else:
                        df = pd.read_pickle(temp_file)
                    
                    # Remove index column if it exists
                    if "index" in df.columns:
                        df = df.drop('index', axis=1)
                    
                    # Create temporary ROOT file
                    import tempfile
                    root_temp = tempfile.NamedTemporaryFile(suffix='.root', delete=False)
                    root_temp.close()
                    
                    with uproot.recreate(root_temp.name) as root_file:
                        root_file[self._region] = df
                    
                    batch_root_files.append(root_temp.name)
                    del df
                    gc.collect()
                    
                except Exception as e:
                    print(f"XGB WARNING: Failed to process temp file {temp_file}: {e}")
                    continue
            
            # Combine this batch into a single ROOT file
            if batch_root_files:
                batch_output = tempfile.NamedTemporaryFile(suffix='.root', delete=False)
                batch_output.close()
                
                hadd_command = f"hadd -f {batch_output.name} " + " ".join(batch_root_files)
                
                import subprocess
                if HADD_AVAILABLE:
                    try:
                        result = subprocess.run(hadd_command, shell=True, capture_output=True, text=True)
                        if result.returncode == 0:
                            root_temp_files.append(batch_output.name)
                            print(f"XGB INFO: Combined batch {batch_start//batch_size + 1}/{(len(temp_files) + batch_size - 1)//batch_size}")
                        else:
                            print(f"XGB ERROR: hadd failed for batch: {result.stderr}")
                            # Fallback: combine manually for this batch
                            self._manual_combine_batch(batch_root_files, batch_output.name)
                            root_temp_files.append(batch_output.name)
                    except Exception as e:
                        print(f"XGB ERROR: Failed to run hadd for batch: {e}")
                        # Fallback: combine manually for this batch
                        self._manual_combine_batch(batch_root_files, batch_output.name)
                        root_temp_files.append(batch_output.name)
                else:
                    print("XGB WARNING: hadd not available, using manual combination for batch...")
                    self._manual_combine_batch(batch_root_files, batch_output.name)
                    root_temp_files.append(batch_output.name)
                
                # Clean up batch temp files
                for batch_file in batch_root_files:
                    try:
                        os.remove(batch_file)
                    except:
                        pass
        
        # Final combination of all batch files
        if root_temp_files:
            print("XGB INFO: Final combination of all batches...")
            
            if HADD_AVAILABLE:
                final_hadd_command = f"hadd -f {output_path} " + " ".join(root_temp_files)
                
                import subprocess
                try:
                    result = subprocess.run(final_hadd_command, shell=True, capture_output=True, text=True)
                    if result.returncode == 0:
                        print("XGB INFO: Successfully combined all files using hadd")
                    else:
                        print(f"XGB ERROR: Final hadd failed: {result.stderr}")
                        # Fallback to manual combination
                        self._fallback_combine_files(root_temp_files, output_path)
                except Exception as e:
                    print(f"XGB ERROR: Failed to run final hadd: {e}")
                    # Fallback to manual combination
                    self._fallback_combine_files(root_temp_files, output_path)
            else:
                print("XGB WARNING: hadd not available, using manual combination...")
                self._fallback_combine_files(root_temp_files, output_path)
            
            # Clean up batch files
            for root_file in root_temp_files:
                try:
                    os.remove(root_file)
                except:
                    pass
    
    def _manual_combine_batch(self, root_files, output_path):
        """Manually combine a batch of ROOT files."""
        if not root_files:
            return
        
        # Read first file
        with uproot.open(root_files[0]) as first_file:
            combined_data = first_file[self._region].arrays(library='pd')
        
        # Append remaining files
        for root_file in root_files[1:]:
            try:
                with uproot.open(root_file) as file:
                    data = file[self._region].arrays(library='pd')
                    combined_data = pd.concat([combined_data, data], ignore_index=True, sort=False)
                    del data
                    gc.collect()
            except Exception as e:
                print(f"XGB WARNING: Failed to read {root_file}: {e}")
                continue
        
        # Write combined result
        with uproot.recreate(output_path) as output_file:
            output_file[self._region] = combined_data
        
        del combined_data
        gc.collect()

def main():

    args=getArgs()
    
    configPath = args.config
    xgb = ApplyXGBHandler(configPath, args.region)

    xgb.setInputFolder(args.inputFolder)
    xgb.setModelFolder(args.modelFolder)

    xgb.loadModels()
    xgb.loadTransformer()

    with open('data/inputs_config.json') as f:
        config = json.load(f)
    sample_list = config['sample_list']

    xgb.setOutputFolder(args.outputFolder)
    for category in sample_list:
        if args.category and category not in args.category: continue
        xgb.applyBDT(category)
    
    xgb.setOutputFolder(args.outputFolder.replace('test', 'val'))
    for category in sample_list:
        if args.category and category not in args.category: continue
        xgb.applyBDT(category, shift=1)

    # xgb.setOutputFolder(args.outputFolder.replace('test', 'train1'))
    # for category in sample_list:
    #     if args.category and category not in args.category: continue
    #     xgb.applyBDT(category, shift=2)

    # xgb.setOutputFolder(args.outputFolder.replace('test', 'train2'))
    # for category in sample_list:
    #     if args.category and category not in args.category: continue
    #     xgb.applyBDT(category, shift=3)
        
    return

if __name__ == '__main__':
    main()
