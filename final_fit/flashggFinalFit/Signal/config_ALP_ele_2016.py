# Config file: options for signal fitting

_year = '2016'

signalScriptCfg = {
  
  # Setup
  # 'inputWSDir':'/vols/cms/jl2117/hgg/ws/UL/Sept20/MC_final/signal_%s'%_year,
  'inputWSDir':'/publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL/output/sig/ele',
  'procs':'auto', # if auto: inferred automatically from filenames
  'cats':'auto', # if auto: inferred automatically from (0) workspace
  'ext':'dcb_%s'%_year,
  'analysis':'STXS', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'120,125,130',
  'ALPmassPoints':'',

  #Photon shape systematics  
  'scales':'HighR9EB,HighR9EE,LowR9EB,LowR9EE,Gain1EB,Gain6EB', # separate nuisance per year
  'scalesCorr':'MaterialCentralBarrel,MaterialOuterBarrel,MaterialForward,FNUFEE,FNUFEB,ShowerShapeHighR9EE,ShowerShapeHighR9EB,ShowerShapeLowR9EE,ShowerShapeLowR9EB', # correlated across years
  'scalesGlobal':'NonLinearity,Geant4', # affect all processes equally, correlated across years
  'smears':'HighR9EBPhi,HighR9EBRho,HighR9EEPhi,HighR9EERho,LowR9EBPhi,LowR9EBRho,LowR9EEPhi,LowR9EERho', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'hep.q'
  #'batch':'condor', # ['condor','SGE','IC','local']
  #'queue':'espresso',

}
