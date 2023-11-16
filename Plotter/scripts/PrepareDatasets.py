#!/usr/bin/python

import os
import optparse
import time
from time import gmtime, strftime
import math

#define function for parsing options
def parseOptions():
    global observalbesTags, modelTags, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-i', '--input', dest='INPUT', type='string',default='./Data/2017/datasets.txt', help='the path of input datastes file')
    parser.add_option('-o', '--output', dest='OUTPUT', type='string',default='./Data/2017', help='the path of output datastes file')
    parser.add_option('--hadd', dest='HADD', action='store_true', default=False , help='hadd signal?')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quiet=0):
    
    print("process cmd:", cmd)
    if not quiet:
        output = os.popen(cmd)
        print("CMD OUTPUT: ", output.read())

def PrepareDatasets():

    # parse the arguments and options
    global opt, args
    parseOptions()

    File_input = opt.INPUT
    File_output = opt.OUTPUT
    File_hadd = opt.HADD

    in_file = open(File_input)

    samples = {}
    for line in in_file:
        if line[0] == "#": continue

        sample_name = line.split()[0]
        sample_path = line.split()[1]
        samples[sample_name] = sample_path
        print("Start to convert sample:", sample_name)

        cmd = 'python3 -m parquet_to_root {input_file} {output_file} -t passedEvents'.format(input_file=sample_path, output_file="{}/{}.root".format(File_output,sample_name))
        processCmd(cmd, quiet=0)

    if 'WminusH' in samples.keys() and 'WplusH' in samples.keys():
        cmd = 'hadd -f {outpath}/WH.root {input_file}'.format(outpath=File_output, input_file=cmd_input.join(['{}/{}.root'.format(File_output, 'WminusH'), '{}/{}.root'.format(File_output, 'WplusH')]))
        processCmd(cmd, quiet=0)

    if File_hadd:
        print('\n')
        print('Start to hadd root files')
        bkgs = []
        sigs = []
        for sample in samples.keys():
            if 'H' in sample:
                sigs.append("{}/{}.root".format(File_output,sample))
            elif 'data' not in sample:
                bkgs.append("{}/{}.root".format(File_output,sample))
        cmd_input = ' '

        cmd = 'hadd -f {outpath}/sig.root {input_file}'.format(outpath=File_output, input_file=cmd_input.join(sigs))
        processCmd(cmd, quiet=0)
        cmd = 'hadd -f {outpath}/bkg.root {input_file}'.format(outpath=File_output, input_file=cmd_input.join(bkgs))
        processCmd(cmd, quiet=0)


# run the submitAnalyzer() as main()
if __name__ == "__main__":
    PrepareDatasets()