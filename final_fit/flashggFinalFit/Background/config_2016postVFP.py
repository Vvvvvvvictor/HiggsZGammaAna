# Config file: options for signal fitting

_year = '2016postVFP'

backgroundScriptCfg = {
  
  # Setup
  'inputWS':f'/eos/user/j/jiehan/root/ws_data/Data_{_year}/Data_{_year}.root', # Input data workspace
  'cats':'auto', # auto: automatically inferred from input ws
  'catOffset':0, # add offset to category numbers (useful for categories from different allData.root files)  
  'ext':f'{_year}', # extension to add to output directory
  'year':f'{_year}', # Use combined when merging all years in category (for plots)

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'hep.q' # for condor e.g. microcentury
  
}
