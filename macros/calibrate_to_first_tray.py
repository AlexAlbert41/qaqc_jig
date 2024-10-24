#! /usr/bin/env python3

import os
import shutil
import glob
import math
import array
import sys
import time

import ROOT
import tdrstyle

import json
import numpy as np
#from typing import NamedTuple


data_path = '/data/QAQC_SM/qaqc-gui_output/SM_QAQC_Production'
selections = []
plotDir = '/data/QAQC_SM/qaqc-gui_output/Calibration_Facors/'


#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.25,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)



modules = []
params = {}

# retrieving root files 
inputFiles = glob.glob(data_path+'/run*/*_analysis.root')
gui_settings_files = glob.glob(data_path+'/run*/qaqc_gui.settings')
print(inputFiles)
for inputFile in inputFiles:
    tokens = inputFile.split('/')
    run = ''
    for token in tokens:
        if 'module' in token:
            module = token[7:21] # SM ID
        if 'run' in token:
            run = int(token[3:]) # run number
    modules.append(module)
    params[module] = [inputFile,run,'GOOD']
#print(modules)

if not os.path.isdir(plotDir):
    os.mkdir(plotDir)


# selecting the modules to be included in the summary: accept 1 if included, 0 otherwise
#print(modules = {}
slot_data_spe_L = {}
slot_data_spe_R = {}
slot_data_src_L = {}
slot_data_src_R = {}
slot_data_LO = {}


slot_data_spe_L_err = {}
slot_data_spe_R_err = {}
slot_data_src_L_err = {}
slot_data_src_R_err = {}
slot_data_LO_err = {}

modules_tested = ["32110020000016"]
for num, module in enumerate(modules):
    #print("Module: ", module)
    if module in modules_tested:
        continue
    modules_tested.append(module)
    param = params[module]
    accept = 1
    if "32110020008456" in module:
    #    print("skipping module")
        continue
    '''
    if param[1] < 50: # measurements re-done after changing module boards
        accept = 0 
    elif param[1] > 56 and param[1] < 63: # uniformity tests
        accept = 0 
    elif param[1]>67:
        accept = 0
    '''
    slot=-1
    for infile in gui_settings_files:
        with open(infile) as myfile:
            data = json.load(myfile)
            #print("barcodes: ", data["barcodes"])
            if module in data["barcodes"]:
                slot=np.where(np.array(data["barcodes"])==module)[0][0]
                #print("slot: ", slot)
    #for selection in selections:
    #    tempAccept = 0
    #    for param in params:
    #        if selection in param:
    #            tempAccept = 1
    #    accept *= tempAccept
    if accept == 0:
        continue
    if slot not in list(slot_data_spe_L.keys()):
        slot_data_spe_L[slot] = []
        slot_data_spe_R[slot] = []
        slot_data_src_L[slot] = []
        slot_data_src_R[slot] = []
        slot_data_LO[slot] = []
    
        slot_data_spe_L_err[slot] = []
        slot_data_spe_R_err[slot] = []
        slot_data_src_L_err[slot] = []
        slot_data_src_R_err[slot] = []
        slot_data_LO_err[slot] = []
    
    
    rootfile = ROOT.TFile(params[module][0],'READ')
    
    # compute calibrations for a single module
    graph_right = rootfile.Get('g_spe_R_vs_bar'); graph_left = rootfile.Get('g_spe_L_vs_bar')
    module_left_spe, module_right_spe = [],[]
    if graph_left.GetN()!=16 or graph_right.GetN()!=16:
        print("One of the SPE graphs does not have 16 channels! Skipping module {}".format(module))
        continue
    for point in range(graph_left.GetN()):
        if graph_left.GetPointY(point) > 2.8 and graph_left.GetPointY(point)<4: 
            module_left_spe.append(graph_left.GetPointY(point))
        else:
            module_left_spe.append(np.nan)
        if graph_right.GetPointY(point) > 2.8 and graph_right.GetPointY(point)<4: 
            module_right_spe.append(graph_right.GetPointY(point))
        else:
            module_right_spe.append(np.nan)
    module_left_spe_avg = np.nanmean(np.array(module_left_spe))
    slot_data_spe_L[slot].append(module_left_spe_avg/module_left_spe) # add scale factors for just one module to dict
    module_right_spe_avg = np.nanmean(np.array(module_right_spe))
    slot_data_spe_R[slot].append(module_right_spe_avg/module_right_spe) # add scale factors for just one module to dict

    
    graph_right = rootfile.Get('g_lyso_R_pc_per_kev_vs_bar'); graph_left = rootfile.Get('g_lyso_L_pc_per_kev_vs_bar')
    module_left_src, module_right_src = [],[]
    if graph_left.GetN()!=16 or graph_right.GetN()!=16:
        print("One of the SRC graphs does not have 16 channels! Skipping module {}".format(module))
        continue
    for point in range(graph_left.GetN()):
        if graph_left.GetPointY(point) > 1.5 and graph_left.GetPointY(point)<2.5: 
            module_left_src.append(graph_left.GetPointY(point))
        else:
            module_left_src.append(np.nan)
        if graph_right.GetPointY(point) > 1.5 and graph_right.GetPointY(point)<2.5: 
            module_right_src.append(graph_right.GetPointY(point))
        else:
            module_right_src.append(np.nan)
    module_left_src_avg = np.nanmean(np.array(module_left_src))
    slot_data_src_L[slot].append(module_left_src_avg/module_left_src) # add scale factors for just one module to dict
    module_right_src_avg = np.nanmean(np.array(module_right_src))
    slot_data_src_R[slot].append(module_right_src_avg/module_right_src) # add scale factors for just one module to dict

hists_dict = {}
print("modules tested: ", len(modules_tested))

for slot in range(12):
    spe_L = np.nanmean(np.array(slot_data_spe_L[slot]), axis=0)
    spe_R = np.nanmean(np.array(slot_data_spe_R[slot]), axis=0)
    spe_total = np.concatenate([spe_L, spe_R])
    temp_spe_hist = ROOT.TGraph()
    for point in range(len(spe_total)):
        temp_spe_hist.SetPointX(point,point)
        temp_spe_hist.SetPointY(point, spe_total[point])
    c = ROOT.TCanvas("SPE_Calib_Factors_slot"+str(slot), "SPE_Calib_Factors_slot"+str(slot), 600, 500)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    temp_spe_hist.SetTitle(";channel;SPE Calibration Ratio")
    temp_spe_hist.SetMarkerStyle(20)
    temp_spe_hist.SetMarkerSize(1.2)
    temp_spe_hist.SetMarkerColor(ROOT.kRed)
    temp_spe_hist.SetLineColor(ROOT.kRed)
    temp_spe_hist.GetYaxis().SetRangeUser(0.9,1.3)
    temp_spe_hist.Draw("PLA")
    c.Print("{}/{}.png".format(plotDir, c.GetName()))
    hists_dict["g_spe_vs_ch_slot{}".format(str(slot))] = temp_spe_hist.Clone()
    del(c)


for slot in range(12):
    src_L = np.nanmean(np.array(slot_data_src_L[slot]), axis=0)
    src_R = np.nanmean(np.array(slot_data_src_R[slot]), axis=0)
    src_total = np.concatenate([src_L, src_R])
    temp_src_hist = ROOT.TGraph()
    for point in range(len(src_total)):
        temp_src_hist.SetPointX(point,point)
        temp_src_hist.SetPointY(point, src_total[point])
    c = ROOT.TCanvas("SRC_Calib_Factors_slot"+str(slot), "SRC_Calib_Factors_slot"+str(slot), 600, 500)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    temp_src_hist.SetTitle(";channel;SRC Calibration Ratio")
    temp_src_hist.SetMarkerStyle(20)
    temp_src_hist.SetMarkerSize(1.2)
    temp_src_hist.SetMarkerColor(ROOT.kRed)
    temp_src_hist.SetLineColor(ROOT.kRed)
    temp_src_hist.GetYaxis().SetRangeUser(0.9,1.3)
    temp_src_hist.Draw("PLA")
    c.Print("{}/{}.png".format(plotDir, c.GetName()))
    hists_dict["g_lyso_pc_per_kev_vs_ch_slot{}".format(str(slot))] = temp_src_hist.Clone()
    del(c)

# creating outfile
outfile = ROOT.TFile("{}/CERN_Module_calib.root".format(plotDir), "RECREATE"    )
outfile.cd()
for key,value in hists_dict.items():
    value.Write(key)
outfile.Close() 



