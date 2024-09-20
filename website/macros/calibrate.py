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




data_path = '/data1/SMQAQC/PRODUCTION/'
selections = ['GOOD']
plotDir = '/data1/html/data1/SMQAQC/PRODUCTION/calibrationPlots_run0050-run0056_calib/'
#plotDir = '/data1/html/data1/SMQAQC/PRODUCTION/calibrationPlots_run0057-run0062/'

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
        #if graph.GetPointY(point) > 1000. and  graph.GetPointY(point) < 5000.:
        htemp.Fill(graph.GetPointY(point))
    return (htemp.GetMean(),htemp.GetRMS())




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
            run = int(token[3:])
    jsonFileName = data_path+'run%04d/qaqc_gui.settings'%run
    config = json.load(open(jsonFileName))
    slot = 0
    for barcode in config['barcodes']:
        if barcode == module and config['module_available'][slot] == 1:
            break
        else:
            slot += 1
    modules.append((module,run))
    params[(module,run)] = [inputFile,slot,'GOOD']
    #print(module,run,params[(module,run)],config)

bad_modules = []
bad_modules.append('32110020000004')
bad_modules.append('32110020000005')
bad_modules.append('32110020000009')
bad_modules.append('32110020000018')
bad_modules.append('32110020000022')
bad_modules.append('32110020000024')
bad_modules.append('32110020000025')
bad_modules.append('32110020000026')
bad_modules.append('32110020000034')
bad_modules.append('32110020000035')
bad_modules.append('32110020000037')
bad_modules.append('32110020000040')




p_spe_L_vs_slot = ROOT.TProfile('p_spe_L_vs_slot','',6,-0.5,5.5)
p_spe_R_vs_slot = ROOT.TProfile('p_spe_R_vs_slot','',6,-0.5,5.5)
p_spe_vs_ampli = ROOT.TProfile('p_spe_vs_ampli','',16,-0.5,15.5)

p_lyso_L_vs_slot = ROOT.TProfile('p_lyso_L_vs_slot','',6,-0.5,5.5)
p_lyso_R_vs_slot = ROOT.TProfile('p_lyso_R_vs_slot','',6,-0.5,5.5)
p_lyso_vs_ampli = ROOT.TProfile('p_lyso_vs_ampli','',16,-0.5,15.5)


print('***** spe *****')
for key in modules:
    module = key[0]
    run = key[1]
    param = params[key]
    slot = param[1]
    accept = 1
    #for selection in selections:
    #    tempAccept = 0
    #    for param in params[module]:
    #        if selection in param:
    #            tempAccept = 1
    #    accept *= tempAccept
    if run < 50 or run > 56:
        accept = 0
    #if run < 57:
    #    accept = 0
    if accept == 0:
        continue
    print(module,run,param)
    
    rootfile = ROOT.TFile(params[(module,run)][0],'READ')
    
    graph = rootfile.Get('g_spe_L_vs_bar')
    mean = GetMeanRMS(graph)[0]
    p_spe_L_vs_slot.Fill(slot,mean)
    
    graph = rootfile.Get('g_spe_R_vs_bar')
    mean = GetMeanRMS(graph)[0]
    p_spe_R_vs_slot.Fill(slot,mean)
    
    graph = rootfile.Get('g_spe_vs_ch')
    mean = GetMeanRMS(graph)[0]
    for point in range(graph.GetN()):
        ch = graph.GetPointX(point)
        if ch < 16:
            p_spe_vs_ampli.Fill(ch%8,graph.GetPointY(point)/mean)
        else:
            p_spe_vs_ampli.Fill(15-ch%8,graph.GetPointY(point)/mean)


print('***** lyso *****')
for key in modules:
    module = key[0]
    run = key[1]
    param = params[key]
    slot = param[1]
    accept = 1
    #for selection in selections:
    #    tempAccept = 0
    #    for param in params[module]:
    #        if selection in param:
    #            tempAccept = 1
    #    accept *= tempAccept
    if module in bad_modules:
        accept = 0
    if run < 50 or run > 56:
        accept = 0
    #if run < 57:
    #    accept = 0
    if accept == 0:
        continue
    print(module,run,param)
    
    rootfile = ROOT.TFile(params[(module,run)][0],'READ')
    
    graph = rootfile.Get('g_lyso_L_pc_per_kev_vs_bar')
    mean = GetMeanRMS(graph)[0]
    p_lyso_L_vs_slot.Fill(slot,mean)
    
    graph = rootfile.Get('g_lyso_R_pc_per_kev_vs_bar')
    mean = GetMeanRMS(graph)[0]
    p_lyso_R_vs_slot.Fill(slot,mean)
    
    graph = rootfile.Get('g_lyso_pc_per_kev_vs_ch')
    mean = GetMeanRMS(graph)[0]
    for point in range(graph.GetN()):
        ch = graph.GetPointX(point)
        if ch < 16:
            p_lyso_vs_ampli.Fill(ch%8,graph.GetPointY(point)/mean)
        else:
            p_lyso_vs_ampli.Fill(15-ch%8,graph.GetPointY(point)/mean)




c = ROOT.TCanvas('c_spe_vs_slot','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
p_spe_L_vs_slot.Scale(1./p_spe_L_vs_slot.GetBinContent(p_spe_L_vs_slot.FindBin(2.)))
p_spe_R_vs_slot.Scale(1./p_spe_R_vs_slot.GetBinContent(p_spe_R_vs_slot.FindBin(2.)))
p_spe_L_vs_slot.SetTitle(';slot;spe charge [a.u.]')
p_spe_L_vs_slot.GetYaxis().SetRangeUser(0.95,1.05)
p_spe_L_vs_slot.SetMarkerStyle(20)
p_spe_L_vs_slot.SetMarkerSize(1.2)
p_spe_L_vs_slot.SetMarkerColor(ROOT.kRed)
p_spe_L_vs_slot.SetLineColor(ROOT.kRed)
p_spe_L_vs_slot.Draw()
p_spe_R_vs_slot.SetMarkerStyle(20)
p_spe_R_vs_slot.SetMarkerSize(1.2)
p_spe_R_vs_slot.SetMarkerColor(ROOT.kBlue)
p_spe_R_vs_slot.SetLineColor(ROOT.kBlue)
p_spe_R_vs_slot.Draw('same')
c.Print('%s/h_spe_LR_slot.png'%plotDir)


c = ROOT.TCanvas('c_spe_vs_ampli','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
p_spe_vs_ampli.SetTitle(';amplifier channel;a.u.')
p_spe_vs_ampli.GetYaxis().SetRangeUser(0.95,1.05)
p_spe_vs_ampli.Draw()
c.Print('%s/h_spe_LR_ch.png'%plotDir)




c = ROOT.TCanvas('c_lyso_vs_slot','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
p_lyso_L_vs_slot.Scale(1./p_lyso_L_vs_slot.GetBinContent(p_lyso_L_vs_slot.FindBin(2.)))
p_lyso_R_vs_slot.Scale(1./p_lyso_R_vs_slot.GetBinContent(p_lyso_R_vs_slot.FindBin(2.)))
p_lyso_L_vs_slot.SetTitle(';slot;lyso charge [a.u.]')
p_lyso_L_vs_slot.GetYaxis().SetRangeUser(0.90,1.05)
p_lyso_L_vs_slot.SetMarkerStyle(20)
p_lyso_L_vs_slot.SetMarkerSize(1.2)
p_lyso_L_vs_slot.SetMarkerColor(ROOT.kRed)
p_lyso_L_vs_slot.SetLineColor(ROOT.kRed)
p_lyso_L_vs_slot.Draw()
p_lyso_R_vs_slot.SetMarkerStyle(20)
p_lyso_R_vs_slot.SetMarkerSize(1.2)
p_lyso_R_vs_slot.SetMarkerColor(ROOT.kBlue)
p_lyso_R_vs_slot.SetLineColor(ROOT.kBlue)
p_lyso_R_vs_slot.Draw('same')
c.Print('%s/h_lyso_LR_slot.png'%plotDir)


c = ROOT.TCanvas('c_lyso_vs_ampli','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
p_lyso_vs_ampli.SetTitle(';amplifier channel;a.u.')
p_lyso_vs_ampli.GetYaxis().SetRangeUser(0.95,1.05)
p_lyso_vs_ampli.Draw()
c.Print('%s/h_lyso_LR_ch.png'%plotDir)



outfile = ROOT.TFile('%s/calibrate.root'%plotDir,'RECREATE')
outfile.cd()
p_spe_L_vs_slot.Write()
p_spe_R_vs_slot.Write()
p_spe_vs_ampli.Write()
p_lyso_L_vs_slot.Write()
p_lyso_R_vs_slot.Write()
p_lyso_vs_ampli.Write()
outfile.Close()
