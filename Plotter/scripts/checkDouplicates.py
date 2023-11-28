####################################################
####################################################

import os
import sys
import numpy as np

from ROOT import *

import argparse
parser = argparse.ArgumentParser(description="Plotter for HZGamma analysis")
parser.add_argument("-i", "--input", dest="input", default="/afs/cern.ch/work/z/zewang/private/HZGamma/data_run2_2017_nanov9/data.root", help="input dataset")
args = parser.parse_args()

treename = "passedEvents"
oldfile = TFile(args.input,'r')
oldchain = oldfile.Get(treename)

entries = oldchain.GetEntries()
print("Number of Entries:", entries)

event_ids = []
n_removed = 0
for ientry in range(entries):
    nb = oldchain.GetEntry(ientry)

    if (ientry % 10000 == 1):
        print("looking at event %d" %ientry)

    #if ientry > 10: break

    event_id = "{}:{}:{}".format(oldchain.run,oldchain.event,oldchain.luminosityBlock)
    #print(event_id)
    if event_id not in event_ids:
        event_ids.append(event_id)
    else:
        n_removed = n_removed + 1

print("Number of duplicated events: ", n_removed)
