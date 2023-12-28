#!/usr/bin/python

import os
import optparse
import json

#define function for parsing options
def parseOptions():
    global observalbesTags, modelTags, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-i', '--input', dest='INPUT', type='string',default='/afs/cern.ch/work/z/zewang/private/HZGamma/data_run2_2017_nanov2_CR', help='the path of input datastes file')
    parser.add_option('-o', '--output', dest='OUTPUT', type='string',default='./Data/2017', help='the path of output datastes file')
    parser.add_option('--hadd', dest='HADD', action='store_true', default=False , help='hadd signal?')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quiet=0):
    
    #print("process cmd:", cmd)
    if not quiet:
        output = os.popen(cmd)
        #print("CMD OUTPUT: ", output.read())
        return output.read()

def getProcessedLumi(file):
    f = open(file)

    lumi = ""

    for line in f.readlines():
        if "[[INFO]] Processed Lumi:" in line:
            lumi = line.split("Lumi: ")[-1]
    return lumi

def CheckProcessedLumi():
    global opt, args
    parseOptions()
    
    cmd = "ls -d {}/*/".format(opt.INPUT)
    out = processCmd(cmd)

    datasets = {}
    for dataset in out.split():
        datasets_name = dataset.split('/')[-2]
        Jobs = processCmd("ls -d {}/{}/*/".format(opt.INPUT, datasets_name))
        lumis = {}
        for iJob in Jobs.split():
            #if ("job_1/" not in iJob) and ("job_6/" not in iJob):
            #    continue
                

            iJob_log = processCmd("ls {}*.log".format(iJob)).rstrip('\n')
            lumi = json.loads(getProcessedLumi("{}".format(iJob_log)).rstrip('\n').replace("'", '"'))

            #print(list(lumi.keys()), list(lumis.keys()))
            
            for run_new in lumi.keys():
                dict_tmp = {}
                #print(run_new, lumis.keys(), lumi[run_new])
                if run_new in lumis.keys():
                    lumis[run_new] = lumis[run_new] + lumi[run_new]
                else:
                    dict_tmp[run_new] = lumi[run_new]
                    lumis.update(dict_tmp)

        #print(lumis, out)
        with open("{}processedLumi.json".format(out.rstrip('\n')),'w') as file:
            file.write(str(lumis).replace("'", '"'))
        file.close()

CheckProcessedLumi()