import os, glob, sys
from optparse import OptionParser
from models import models

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG T2W RUN II ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")

def get_options():
  parser = OptionParser()

  parser.add_option('--mass_ALP', dest='mass_ALP', default=1, type='int', help="ALP mass") # PZ
  parser.add_option('--year', dest='year', default='2022preEE,2022postEE', help="Comma separated list of years")
  parser.add_option("--channel", dest='channel', default='', help="ele or mu") # PZ

  parser.add_option('--mode', dest='mode', default='mu_inclusive', help="Physics Model (specified in models.py)")
  parser.add_option('--ext',dest='ext', default="", help='In case running over datacard with extension')
  parser.add_option('--common_opts',dest='common_opts', default="-m 125 higgsMassRange=122,128", help='Common options')
  parser.add_option('--batch', dest='batch', default='local', help="Batch system [SGE,IC,condor]")
  parser.add_option('--queue', dest='queue', default='workday', help="Condor queue")
  parser.add_option('--ncpus', dest='ncpus', default=4, type='int', help="Number of cpus")
  parser.add_option('--dryRun', dest='dryRun', action="store_true", default=True, help="Only create submission files")
  return parser.parse_args()
(opt,args) = get_options()

def leave():
  print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG T2W RUN II (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit(1)

def run(cmd):
  print("%s\n\n"%cmd)
  os.system(cmd)

if opt.mode not in models: 
  print(" --> [ERROR] opt.mode (%s) is not specified in models.py. Leaving..."%opt.mode)
  leave()

extStr = "_%s"%opt.channel if opt.channel != '' else ''
print(" --> Running text2workspace for model: %s"%opt.mode)
print(" --> Input: /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard%s/%s_pruned_datacard_%s_%s.txt --> Output: ./output_datacard_rootfile%s/%s_Datacard_%s_%s_%s.root"%(extStr,opt.mass_ALP,opt.year,opt.channel, extStr,opt.mass_ALP,opt.year,opt.channel,opt.mode,))

if not os.path.isdir(f"./t2w_jobs_{opt.channel}"): os.system(f"mkdir ./t2w_jobs_{opt.channel}")
if not os.path.isdir(f"./output_datacard_rootfile_{opt.channel}"): os.system(f"mkdir ./output_datacard_rootfile_{opt.channel}")

# Open submission file to write to
fsub = open("./t2w_jobs_%s/%s_t2w_%s_%s_%s.sh"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel),"w")
fsub.write("#!/bin/bash\n\n")
fsub.write("cd %s\n\n"%os.environ['PWD'])
fsub.write("eval `scramv1 runtime -sh`\n\n")
# if not os.path.isdir("./output_Datacard%s"%extStr): os.system("mkdir ./output_Datacard%s"%extStr)
# fdataName = "./output_Datacard%s/%s_pruned_datacard_%s_%s.txt"%(extStr,opt.mass_ALP,opt.year,opt.channel)
# fsub.write("text2workspace.py Datacard%s.txt -o Datacard%s_%s.root %s %s"%(opt.ext,opt.ext,opt.mode,opt.common_opts,models[opt.mode]))
fsub.write("text2workspace.py /publicfs/cms/user/laipeizhu/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard%s/%s_pruned_datacard_%s_%s.txt -o ./output_datacard_rootfile%s/%s_Datacard_%s_%s_%s.root %s %s"%(extStr,opt.mass_ALP,opt.year,opt.channel,  extStr,opt.mass_ALP,opt.year,opt.channel,opt.mode,  opt.common_opts,models[opt.mode]))
fsub.close()

# Change permission for file
os.system("chmod 775 ./t2w_jobs_%s/%s_t2w_%s_%s_%s.sh"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel))

# If using condor then also write submission file
if opt.batch == 'condor':
  f_cdr = open("./t2w_jobs_%s/%s_t2w_%s_%s_%s.sub"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel),"w")
  f_cdr.write("executable          = %s/src/flashggFinalFit/Combine/t2w_jobs_%s/%s_t2w_%s_%s_%s.sh\n"%(os.environ['CMSSW_BASE'],opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel))
  f_cdr.write("output              = %s/src/flashggFinalFit/Combine/t2w_jobs_%s/%s_t2w_%s_%s_%s.sh.out\n"%(os.environ['CMSSW_BASE'],opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel))
  f_cdr.write("error               = %s/src/flashggFinalFit/Combine/t2w_jobs_%s/%s_t2w_%s_%s_%s.sh.err\n"%(os.environ['CMSSW_BASE'],opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel))
  f_cdr.write("log                 = %s/src/flashggFinalFit/Combine/t2w_jobs_%s/%s_t2w_%s_%s_%s.sh.log\n"%(os.environ['CMSSW_BASE'],opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel))
  f_cdr.write("+JobFlavour         = \"%s\"\n"%opt.queue)
  f_cdr.write("RequestCpus         = %g\n"%opt.ncpus)
  f_cdr.write("queue\n")
  f_cdr.close()

# Submit
if opt.batch == "condor": 
  if os.environ['PWD'].startswith("/eos"):
    subcmd = "condor_submit -spool ./t2w_jobs_%s/%s_t2w_%s_%s_%s.sub"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel) 
  else:
    subcmd = "condor_submit ./t2w_jobs_%s/%s_t2w_%s_%s_%s.sub"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel)
elif opt.batch == 'local': subcmd = "bash ./t2w_jobs_%s/%s_t2w_%s_%s_%s.sh"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel)
else: subcmd = "qsub -q hep.q -l h_rt=6:0:0 -l h_vmem=24G ./t2w_jobs_%s/%s_t2w_%s_%s_%s.sh"%(opt.channel,opt.mass_ALP,opt.mode,opt.year,opt.channel)
if opt.dryRun: print("[DRY RUN] %s"%subcmd)
else: run(subcmd)
