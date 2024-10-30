# Config file: options for signal fitting

_year = '2017'

signalScriptCfg = {
  
  # Setup
  # 'inputWSDir':'PATH_TO_INPUTS/workspaces/signal_%s'%_year,
  'inputWSDir':'/eos/user/j/jiehan/root/ws_cor_syst/ggH/2017/',
  'procs':'auto', # if auto: inferred automatically from filenames
  'cats':'auto', # if auto: inferred automatically from (0) workspace
  'ext':'tutorial_%s'%_year,
  'analysis':'tutorial', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'120,125,130',

  #Photon shape systematics  
  'scales':'scale', # separate nuisance per year
  'scalesCorr':'', # correlated across years
  'scalesGlobal':'', # affect all processes equally, correlated across years
  'smears':'smear', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'espresso',

}
