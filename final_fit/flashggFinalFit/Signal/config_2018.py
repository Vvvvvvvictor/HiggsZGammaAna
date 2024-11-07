# Config file: options for signal fitting

_year = '2018'

signalScriptCfg = {
  
  # Setup
  'inputWSDir':'/eos/user/j/jiehan/root/ws_cor_syst/signal_%s/'%_year, #i.e. 'PATH_TO_INPUTS/workspaces/signal_%s'%_year
  'procs':'auto', # if auto: inferred automatically from filenames
  'cats':'auto', # if auto: inferred automatically from (0) workspace
  'ext':'%s'%_year,
  'analysis':'tutorial', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'125', #'120,125,130'

  #Photon shape systematics  
  'scales':'scale', # separate nuisance per year
  'scalesCorr':'fnuf,material', # correlated across years
  'scalesGlobal':'', # affect all processes equally, correlated across years
  'smears':'smear', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'espresso',

}
