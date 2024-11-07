# Input config file for running trees2ws

trees2wsCfg = {

  # Name of RooDirectory storing input tree
  'inputTreeDir':'',

  # Variables to be added to dataframe: use wildcard * for common strings
  'mainVars':["CMS_hzg_mass","weight", "weight_*"], # Vars to add to nominal RooDatasets ‘dZ’ may be necessary
  'dataVars':["CMS_hzg_mass","weight"], # Vars for data workspace (trees2ws_data.py script)
  'stxsVar':'', # Var for STXS splitting: if using option doSTXSSplitting
  'systematicsVars':[], # Variables to add to sytematic RooDataHists
  'theoryWeightContainers':{}, # Theory weights to add to nominal + NOTAG RooDatasets, value corresponds to number of weights (0-N)

  # List of systematics: use string YEAR for year-dependent systematics
  'systematics':[],

  # Analysis categories: python list of cats or use 'auto' to extract from input tree
  'cats': ["ggH0", "ggH1", "ggH2", "ggH3", "VBF0", "VBF1", "VBF2", "VBF3", "lep", "VH", "ZH", "ttHh", "ttHl"]

}
