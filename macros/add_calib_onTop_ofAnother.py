#! /usr/bin/env python3
import os
import shutil
import glob
import math
import array
import sys
import time
import json

import ROOT
import tdrstyle


# paths 
calib_path = '/home/cmsdaq/DAQ/qaqc_calibration/'
calib1 = '{}/master_calib.root'.format(calib_path)
calib2 = '{}/digitizers_calib.root'.format(calib_path)
out_calib = '{}/master_calib_withDigi.root'.format(calib_path)


# graphs 
graphs1 = ['g_spe_vs_ch', "g_lyso_pc_per_kev_vs_ch"]
graphs2 = ['g_spe', "g_lyso"]

spe = {}
lyso = {}

spe2 = ROOT.TFile(calib2,"OPEN").Get("g_spe_vs_ch_average")
lyso2 = ROOT.TFile(calib2,"OPEN").Get("g_lyso_pc_per_kev_vs_ch_average") 

for slot in range(6):
    spe1 = ROOT.TFile(calib1,"OPEN").Get("g_spe_slot{}".format(slot))
    lyso1 = ROOT.TFile(calib1,"OPEN").Get("g_lyso_slot{}".format(slot))

    spe[slot] = ROOT.TGraphErrors()
    lyso[slot] = ROOT.TGraphErrors()

    for ch in range(spe1.GetN()):
        val1 = spe1.GetY()[ch]
        val2 = spe2.GetY()[ch]
        spe[slot].SetPoint(spe[slot].GetN(), ch, val1*val2)

        val1 = lyso1.GetY()[ch]
        val2 = lyso2.GetY()[ch]
        lyso[slot].SetPoint(lyso[slot].GetN(), ch, val1*val2)

# creating outfile
outfile = ROOT.TFile("{}/master_calib_with_digitizers.root".format(calib_path), "RECREATE")
outfile.cd()
for slot in range(6):
    spe[slot].Write("g_spe_slot{}".format(slot))
    lyso[slot].Write("g_lyso_slot{}".format(slot))
outfile.Close()
