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
data_path = '/data1/SMQAQC/PRODUCTION/'
calib_path = '/home/cmsdaq/DAQ/qaqc_calibration/'
plotDir = '/data1/html/data1/SMQAQC/PRODUCTION/calibrationPlots_digitizers_run0079-run0083/'

if not os.path.isdir(plotDir):
    os.mkdir(plotDir)

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.25,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)

def GetMeanRMS(graph):
    htemp = ROOT.TH1F('htemp','',100,0.,10000)
    for point in range(graph.GetN()):
        htemp.Fill(graph.GetPointY(point))
    return (htemp.GetMean(),htemp.GetRMS())


digis = {
    "crate": data_path+'/run0079/' ,
    "desktop" : data_path+'/run0083/' }

modules = {
    "70": 'module_32110020000070_analysis.root',
    "71": 'module_32110020000071_analysis.root',
    "72": 'module_32110020000072_analysis.root',
}

graph_names = {
    'g_spe_vs_ch': "spe charge [a.u.]",
    'g_lyso_pc_per_kev_vs_ch': "lyso charge [a.u.]"}

# retrieve graphs from root files
graphs = {}
for g_name in graph_names:
    graphs[g_name] = {}
    for module in modules:
        graphs[g_name][module] = {}
        for digi in digis:
            graphs[g_name][module][digi] = ROOT.TFile("{}/{}".format(digis[digi],modules[module]),"OPEN").Get(g_name)

# create graphs ratio
graphs_ratio = {}
graph_ave = {}
for g_name in graph_names:
    graphs_ratio[g_name] = {}
    graph_ave[g_name] = ROOT.TGraphErrors()
    for module in modules:
        graphs_ratio[g_name][module] = ROOT.TGraphErrors()

# fill graphs
for g_name in graph_names:
    for module in modules:
        channels = graphs[g_name][module]["crate"].GetN()
        for ch in range(channels):
            ratio = graphs[g_name][module]["crate"].GetY()[ch] / graphs[g_name][module]["desktop"].GetY()[ch]
            graphs_ratio[g_name][module].SetPoint(graphs_ratio[g_name][module].GetN(), ch, ratio)
            
# compute average between 3 modules
for g_name in graph_names:
    channels = graphs[g_name]["70"]["crate"].GetN() 
    for ch in range(channels): 
        mean = graphs_ratio[g_name]["70"].GetY()[ch] + graphs_ratio[g_name]["71"].GetY()[ch] + graphs_ratio[g_name]["72"].GetY()[ch]
        mean /= 3
        graph_ave[g_name].SetPoint(graph_ave[g_name].GetN(), ch, mean)

# draw
for g_name in graph_names:
    for module in modules:
        c = ROOT.TCanvas(g_name+"_"+module, g_name+"_"+module, 600, 500)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        graphs_ratio[g_name][module].SetTitle(";channel;{}".format(graph_names[g_name]))
        graphs_ratio[g_name][module].SetMarkerStyle(20)
        graphs_ratio[g_name][module].SetMarkerSize(1.2)
        graphs_ratio[g_name][module].SetMarkerColor(ROOT.kRed)
        graphs_ratio[g_name][module].SetLineColor(ROOT.kRed)
        graphs_ratio[g_name][module].GetYaxis().SetRangeUser(0.9,1.3)  
        graphs_ratio[g_name][module].Draw("PLA")
        c.Print("{}/{}.png".format(plotDir, c.GetName()))
        del(c)

    c = ROOT.TCanvas(g_name+"_average", g_name+"_average", 600, 500)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    graph_ave[g_name].SetTitle(";channel;avg. {}".format(graph_names[g_name]))
    graph_ave[g_name].SetMarkerStyle(20)
    graph_ave[g_name].SetMarkerSize(1.2)
    graph_ave[g_name].SetMarkerColor(ROOT.kMagenta)
    graph_ave[g_name].SetLineColor(ROOT.kMagenta)
    graph_ave[g_name].GetYaxis().SetRangeUser(0.9,1.3)  
    graph_ave[g_name].Draw("PLA")
    c.Print("{}/{}.png".format(plotDir, c.GetName()))
    del(c)

# creating outfile
outfile = ROOT.TFile("{}/digitizers_calib.root".format(calib_path), "RECREATE")
outfile.cd()
for g_name in graph_names:
    graph_ave[g_name].Write("{}_average".format(g_name))
outfile.Close()
