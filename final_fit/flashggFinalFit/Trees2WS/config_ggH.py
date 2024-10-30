# Input config file for running trees2ws

trees2wsCfg = {

  # Name of RooDirectory storing input tree
  'inputTreeDir':'',

  # Variables to be added to dataframe: use wildcard * for common strings
  'mainVars':["H_mass","weight", "weight_*"], # Vars to add to nominal RooDatasets ‘dZ’ may be necessary
  'dataVars':["H_mass","weight"], # Vars for data workspace (trees2ws_data.py script)
  'stxsVar':'', # Var for STXS splitting: if using option doSTXSSplitting
  'systematicsVars':["CMS_hgg_mass","weight"], # Variables to add to sytematic RooDataHists
  'theoryWeightContainers':{}, # Theory weights to add to nominal + NOTAG RooDatasets, value corresponds to number of weights (0-N)

  # List of systematics: use string YEAR for year-dependent systematics
  'systematics':["Scale","Smearing"],

  # Analysis categories: python list of cats or use 'auto' to extract from input tree
  'cats': ["zero_to_one_jet"]

}
