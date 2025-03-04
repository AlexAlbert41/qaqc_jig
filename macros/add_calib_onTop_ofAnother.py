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
calib_path = '/data/QAQC_SM/qaqc-gui_output/calibs_sodium'
calib1 = '{}/calib_scale.root'.format(calib_path)
calib2 = '{}/calib_channels.root'.format(calib_path)
out_calib = '{}/master_calib_channels_scale.root'.format(calib_path)


# graphs 
graphs1 = ['g_spe_vs_ch', "g_lyso_pc_per_kev_vs_ch"]
graphs2 = ['g_spe', "g_lyso"]

spe = {}
lyso = {}

#spe2 = ROOT.TFile(calib2,"OPEN").Get("g_spe_vs_ch_average")
#lyso2 = ROOT.TFile(calib2,"OPEN").Get("g_lyso_pc_per_kev_vs_ch_average") 

rootfile = ROOT.TFile(calib1,'READ')
speL = rootfile.Get("p_spe_L_vs_slot")
lysoL = rootfile.Get("p_lyso_L_vs_slot")
speR = rootfile.Get("p_spe_R_vs_slot")
lysoR = rootfile.Get("p_lyso_R_vs_slot")

speL_graph = ROOT.TGraphErrors()
speR_graph = ROOT.TGraphErrors()
lysoL_graph = ROOT.TGraphErrors()
lysoR_graph = ROOT.TGraphErrors()
for i in range(12):
    speL_graph.SetPointX(i, i)
    speL_graph.SetPointY(i, speL.GetBinContent(speL.FindBin(i)))

    speR_graph.SetPointX(i, i)
    speR_graph.SetPointY(i, speR.GetBinContent(speR.FindBin(i)))

    lysoL_graph.SetPointX(i, i)
    lysoL_graph.SetPointY(i, lysoL.GetBinContent(lysoL.FindBin(i)))

    lysoR_graph.SetPointX(i, i)
    lysoR_graph.SetPointY(i, lysoR.GetBinContent(lysoR.FindBin(i)))

for slot in range(12):
    
    spe2 = ROOT.TFile(calib2,"OPEN").Get("g_spe_vs_ch_slot{}".format(slot))
    lyso2 = ROOT.TFile(calib2, "OPEN").Get("g_lyso_pc_per_kev_vs_ch_slot{}".format(slot))
    
    spe[slot] = ROOT.TGraphErrors()
    lyso[slot] = ROOT.TGraphErrors()

    for ch in range(spe2.GetN()):
        val2_spe = spe2.GetY()[ch]
        val2_lyso = lyso2.GetY()[ch]
        
        if ch<16: #left side
            val1_spe = 1/speL_graph.GetY()[slot]
            val1_lyso = 1/lysoL_graph.GetY()[slot]

        else: #right side
            val1_spe = 1/speR_graph.GetY()[slot]
            val1_lyso = 1/lysoR_graph.GetY()[slot]

        spe[slot].SetPoint(spe[slot].GetN(), ch, val1_spe*val2_spe)
        lyso[slot].SetPoint(lyso[slot].GetN(), ch, val1_lyso*val2_lyso)

# creating outfile
outfile = ROOT.TFile("{}/master_calib.root".format(calib_path), "RECREATE")
outfile.cd()
for slot in range(12):
    spe[slot].Write("g_spe_slot{}".format(slot))
    lyso[slot].Write("g_lyso_slot{}".format(slot))
outfile.Close()
