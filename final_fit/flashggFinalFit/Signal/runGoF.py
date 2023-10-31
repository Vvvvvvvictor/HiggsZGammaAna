#!/usr/bin/env python

from optparse import OptionParser
import os
import sys
import copy
import math

import ROOT as r

parser = OptionParser()
parser.add_option("-i",help="Input data file path")
(options, args) = parser.parse_args()

method = ["saturated", "KS"]

for m in method:
    file_obs_name = "higgsCombine_"+m+".GoodnessOfFit.mH125.root"
    file_toy_name = "higgsCombine_"+m+".GoodnessOfFit.mH125.12345.root"

    file_Obs = r.TFile.Open(file_obs_name)
    file_Toy = r.TFile.Open(file_toy_name)

    tree_Obs = file_Obs.Get("limit")
    tree_Toy = file_Toy.Get("limit")

    n_obs = tree_Obs.GetEntriesFast()
    for jentry in range(n_obs):
        nb = tree_Obs.GetEntry(jentry)
        obs = tree_Obs.limit

    n_toy = tree_Toy.GetEntriesFast()
    n_toy_pass = 0
    for ientry in range(n_toy):
        na = tree_Toy.GetEntry(ientry)
        if tree_Toy.limit > obs:
            n_toy_pass = n_toy_pass + 1

    p_value = float(n_toy_pass)/float(n_toy)

    print "Method:", m, ", number of toys:", n_toy, "p_obs:", obs, "p_value:", p_value
