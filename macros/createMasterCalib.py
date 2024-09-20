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

from typing import NamedTuple




inputfile_vs_slot = ROOT.TFile('/data1/html/data1/SMQAQC/PRODUCTION/calibrationPlots_run0057-run0062/calibrate.root','READ')
inputfile_vs_ampli = ROOT.TFile('/data1/html/data1/SMQAQC/PRODUCTION/calibrationPlots_run0050-run0056/calibrate.root','READ')

p_spe_L_vs_slot = inputfile_vs_slot.Get('p_spe_L_vs_slot')
p_spe_R_vs_slot = inputfile_vs_slot.Get('p_spe_R_vs_slot')
p_lyso_L_vs_slot = inputfile_vs_slot.Get('p_lyso_L_vs_slot')
p_lyso_R_vs_slot = inputfile_vs_slot.Get('p_lyso_R_vs_slot')

p_spe_vs_ampli = inputfile_vs_ampli.Get('p_spe_vs_ampli')
p_lyso_vs_ampli = inputfile_vs_ampli.Get('p_lyso_vs_ampli')




outfile = ROOT.TFile('master_calib.root','RECREATE')

g_spe = {}
g_lyso = {}

for slot in range(6):
    
    g_spe[slot] = ROOT.TGraph()
    g_lyso[slot] = ROOT.TGraph()
    
    slot_L_spe_calib = 1. / p_spe_L_vs_slot.GetBinContent(slot+1)
    slot_R_spe_calib = 1. / p_spe_R_vs_slot.GetBinContent(slot+1)
    slot_L_lyso_calib = 1. / p_lyso_L_vs_slot.GetBinContent(slot+1)
    slot_R_lyso_calib = 1. / p_lyso_R_vs_slot.GetBinContent(slot+1)    
    
    for ch in range(32):
        
        if ch < 16:
            ampli_spe_calib = 1. / p_spe_vs_ampli.GetBinContent(ch%8+1)
            ampli_lyso_calib = 1. / p_lyso_vs_ampli.GetBinContent(ch%8+1)
            g_spe[slot].SetPoint(ch,ch,slot_L_spe_calib*ampli_spe_calib)
            g_lyso[slot].SetPoint(ch,ch,slot_L_lyso_calib*ampli_lyso_calib)
        else:
            ampli_spe_calib = 1. / p_spe_vs_ampli.GetBinContent(15-ch%8+1)
            ampli_lyso_calib = 1. / p_lyso_vs_ampli.GetBinContent(15-ch%8+1)
            g_spe[slot].SetPoint(ch,ch,slot_R_spe_calib*ampli_spe_calib)
            g_lyso[slot].SetPoint(ch,ch,slot_R_lyso_calib*ampli_lyso_calib)

outfile.cd()
for slot in range(6):
    g_spe[slot].Write('g_spe_slot%d'%slot)
    g_lyso[slot].Write('g_lyso_slot%d'%slot)
outfile.Close()
