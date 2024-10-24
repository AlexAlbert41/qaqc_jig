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
CERN_module_file = '/home/cptlab/qaqc-gui_output/CERN_Module_CERN_Results/module_32110020000016_analysis.root'
data_path = '/home/cptlab/qaqc-gui_output/modules_newRuns/calibrationData_CERNModule/12_Calibration_Runs/'
calib_output_path = '/home/cptlab/qaqc-gui_output/CERN_Module_Calibration_Output/'
plotDir = '/home/cptlab/qaqc-gui_output/CERN_Module_Calibration_Output'

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




graph_names = {
    'g_spe_vs_ch': "spe charge [a.u.]",
    'g_lyso_pc_per_kev_vs_ch': "lyso charge [a.u.]"}


modules = []
params = {}
inputFiles = glob.glob(data_path+'/run*/*_analysis.root')
for inputFile in inputFiles:
    tokens = inputFile.split('/')
    run = ''
    for token in tokens:
        if 'module' in token:
            module = token[7:21]
        if 'run' in token:
            run = token[3:]
    jsonFileName = data_path+'run%s/qaqc_gui.settings'%run
    config = json.load(open(jsonFileName))
    slot = 0
    for barcode in config['barcodes']:
        if barcode == module and config['module_available'][slot] == 1:
            break
        else:
            slot += 1
    modules.append((module,run))
    params[(module,run)] = [inputFile,slot,'GOOD']


#get graphs from CERN module to calibrate to
calib_file = ROOT.TFile(CERN_module_file, "OPEN")
calib_graphs = {}
for g_name in graph_names:
    calib_graphs[g_name] = calib_file.Get(g_name)


graphs = {}
graphs_ratio = {}
for g_name in graph_names:
    graphs[g_name] = {}
    graphs_ratio[g_name] = {}

for module in modules:
    # retrieve graphs from root files
    slot = params[module][1]
    for g_name in graph_names:
        graphs[g_name][slot] = ROOT.TFile(data_path+"run{}/module_{}_analysis.root".format(module[1],module[0]),"OPEN").Get(g_name)
    # create graphs ratio
    for g_name in graph_names:
        graphs_ratio[g_name][slot] = ROOT.TGraphErrors()

    # fill graphs
    for g_name in graph_names:
        channels = graphs[g_name][slot].GetN()
        for ch in range(channels):
            ratio = calib_graphs[g_name].GetY()[ch] / graphs[g_name][slot].GetY()[ch]
            graphs_ratio[g_name][slot].SetPoint(graphs_ratio[g_name][slot].GetN(), ch, ratio)
                
    # draw
    for g_name in graph_names:
        c = ROOT.TCanvas(g_name+"_slot"+str(slot), g_name+"_slot"+str(slot), 600, 500)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        graphs_ratio[g_name][slot].SetTitle(";channel;{}".format(graph_names[g_name]))
        graphs_ratio[g_name][slot].SetMarkerStyle(20)
        graphs_ratio[g_name][slot].SetMarkerSize(1.2)
        graphs_ratio[g_name][slot].SetMarkerColor(ROOT.kRed)
        graphs_ratio[g_name][slot].SetLineColor(ROOT.kRed)
        graphs_ratio[g_name][slot].GetYaxis().SetRangeUser(0.9,1.3)  
        graphs_ratio[g_name][slot].Draw("PLA")
        c.Print("{}/{}.png".format(plotDir, c.GetName()))
        del(c)

# creating outfile
outfile = ROOT.TFile("{}/CERN_Module_calib.root".format(calib_output_path), "RECREATE")
outfile.cd()
for g_name in graph_names:
    for slot in range(12):
        graphs_ratio[g_name][slot].Write("{}_slot{}".format(g_name, str(slot)))
outfile.Close()