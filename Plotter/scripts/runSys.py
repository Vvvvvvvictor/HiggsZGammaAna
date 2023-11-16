#!/usr/bin/python
import glob
import sys, os, pwd, commands
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
import subprocess

import datetime
import threading
import multiprocessing

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
    #    print cmd
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
        return "ERROR!!! "+output
    else:
        return output

def execCmd(cmd):
    try:
        print "cmd: %s start running%s" % (cmd,datetime.datetime.now())
        #os.system(cmd)
        subprocess.call(cmd)
        print "cmd: %s end running%s" % (cmd,datetime.datetime.now())
    except Exception, e:
        print '%s\t failed,reason: \r\n%s' % (cmd,e)

def process(cmd,q):
    #subprocess.call(cmd)
    print list(cmd)
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    lines = p.stdout.read()
    data = []
    for i in lines.split('\n'):
        if '\t' in i:
            data.append(i)
    q.put(data)

    #for line in p.stdout.read():
    #    if line.startswith('Printing event yield:'):
    #        print line


def hadd():
    start = time.time()
    verbose = 0

    sysNames = ['norm', 'ShowerShape', 'pho_scale', 'pho_smear', 'lep_scale', 'lep_smear']
    q = multiprocessing.Queue()

    cmds = {}
    outputs = []
    sub_p = {}

    for sysName in sysNames:

        #cmds[sysName] = ['python', 'ALP_BDTSys.py', '-C', '-s'] + [sysName]
        #sub_p[sysName] = multiprocessing.Process(target=process,args=(cmds[sysName],q))
        cmds[sysName] = 'python ALP_BDTSys.py -C -s ' + sysName
        output = processCmd(cmds[sysName])
        data = []
        for i in output.split('\n'):
            if '\t' in i:
                data.append(i)
        outputs.append(data)

    print outputs

    '''
    for sysName in sysNames:
        sub_p[sysName].start()

    for sysName in sysNames:
        sub_p[sysName].join()

    results = [q.get() for j in sub_p]

    print results
    '''
    for i in range(len(outputs[0])):
        line = ''
        for j in range(len(outputs)):
            if j != 0:
                line = line + outputs[j][i].replace(outputs[j][i].split('\t\t')[0],'')
            else:
                line = line + outputs[j][i]
        print line

    end = time.time()
    print str(round(end-start,3))+'s'




# run the submitAnalyzer() as main()
if __name__ == "__main__":
    hadd()
