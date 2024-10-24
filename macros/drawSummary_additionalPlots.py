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


data_path = '/data/QAQC_SM/qaqc-gui_output/Summary_Plots_First37_NoCal/SM_QAQC_Production/'
selections = []
plotDir = '/data/QAQC_SM/qaqc-gui_output/Summary_Plots_NewDirectory_NoCal_AddedPlots/'


#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.25,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)

def GetMaxVar(graph):
    minVal = 999999.
    maxVal = -999999.
    for point in range(graph.GetN()):
        if graph.GetPointY(point) < minVal:
            minVal = graph.GetPointY(point)
        if graph.GetPointY(point) > maxVal:
            maxVal = graph.GetPointY(point)            
    return maxVal-minVal


def GetMeanRMS(graph):
    htemp = ROOT.TH1F('htemp','',100,-100.,10000)
    for point in range(graph.GetN()):
        #if graph.GetPointY(point) > 1000. and  graph.GetPointY(point) < 5000.:
        htemp.Fill(graph.GetPointY(point))
    return (htemp.GetMean(),htemp.GetRMS())

def GetMeanRMS_abs(graph):
    htemp = ROOT.TH1F('htemp','',100,-100.,10000)
    for point in range(graph.GetN()):
        #if graph.GetPointY(point) > 1000. and  graph.GetPointY(point) < 5000.:
        htemp.Fill(abs(graph.GetPointY(point)))
    return (htemp.GetMean(),htemp.GetRMS())


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

# creating histos
h_spe_L_ch = ROOT.TH1F('h_spe_L_ch','',100,3.,5.)
h_spe_R_ch = ROOT.TH1F('h_spe_R_ch','',100,3.,5.)

h_LO_avg_bar = ROOT.TH1F('h_LO_avg_bar','',100,1200.,5200.)
h_LO_L_bar = ROOT.TH1F('h_LO_L_bar','',100,1200.,5200.)
h_LO_R_bar = ROOT.TH1F('h_LO_R_bar','',100,1200.,5200.)
h_LO_asymm_bar = ROOT.TH1F('h_LO_asymm_bar','',100,0,0.2)

h_LO_avg_ch = ROOT.TH1F('h_LO_avg_ch','',100,1200.,5200.)
h_LO_L_ch = ROOT.TH1F('h_LO_L_ch','',100,1200.,5200.)
h_LO_R_ch = ROOT.TH1F('h_LO_R_ch','',100,1200.,5200.)
h_LO_asymm_ch = ROOT.TH1F('h_LO_asymm_ch','',100,-0.4,0.4)

h_LOrms_bar = ROOT.TH1F('h_LOrms_bar','',60,0.,30.)
h_LOrms_ch = ROOT.TH1F('h_LOrms_ch','',60,0.,30.)

h_LOmaxvar_bar = ROOT.TH1F('h_LOmaxvar_bar','',50,0.,100.)
h_LOmaxvar_ch = ROOT.TH1F('h_LOmaxvar_ch','',50,0.,100.)

#added graphs for additional analyses
g_spe_vs_LO_ch_L = ROOT.TGraph()
g_spe_vs_LO_ch_R = ROOT.TGraph()
g_spe_vs_asymm_ch_L = ROOT.TGraph()
g_spe_vs_asymm_ch_R = ROOT.TGraph()
g_LO_vs_asymm_ch_R = ROOT.TGraph()
g_LO_vs_asymm_ch_L = ROOT.TGraph()
g_LO_vs_asymm_bar = ROOT.TGraph()
g_spe_vs_src_ch_L = ROOT.TGraph()
g_spe_vs_src_ch_R = ROOT.TGraph()
g_src_vs_LO_ch_L = ROOT.TGraph()
g_src_vs_LO_ch_R = ROOT.TGraph()
g_spe_asymm_vs_LO_asymm = ROOT.TGraph()
g_src_asymm_vs_LO_asymm = ROOT.TGraph()

g_slot_vs_LO = ROOT.TGraphErrors()
g_slot_vs_spe_L = ROOT.TGraphErrors()
g_slot_vs_spe_R = ROOT.TGraphErrors()
g_slot_vs_src_L = ROOT.TGraphErrors()
g_slot_vs_src_R = ROOT.TGraphErrors()

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
    print("Module: ", module)
    if module in modules_tested:
        continue
    modules_tested.append(module)
    param = params[module]
    accept = 1
    #if "32110020008456" in module:
    #    print("skipping module")
    #    continue
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
                print("slot: ", slot)
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
    
    # filling histos
    graph = rootfile.Get('g_spe_L_vs_bar')
    c = ROOT.TCanvas('c', 'c', 800, 800)
    graph.Draw()
    #c.SaveAs(plotdir+str(params[module][0])+"_graph.png")
    for point in range(graph.GetN()):
        h_spe_L_ch.Fill(graph.GetPointY(point))
        g_spe_vs_LO_ch_L.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_asymm_ch_L.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_src_ch_L.SetPointX(num*16+point, graph.GetPointY(point))


    graph = rootfile.Get('g_spe_R_vs_bar'); graph2 = rootfile.Get('g_spe_L_vs_bar')
    for point in range(graph.GetN()):
        h_spe_R_ch.Fill(graph.GetPointY(point))
        g_spe_vs_LO_ch_R.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_asymm_ch_R.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_src_ch_R.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_asymm_vs_LO_asymm.SetPointX(num*16+point, 2*(graph2.GetPointY(point)-graph.GetPointY(point))/(graph2.GetPointY(point)+graph.GetPointY(point)))
    slot_data_spe_L[slot].append(GetMeanRMS(graph2)[0]) 
    slot_data_spe_R[slot].append(GetMeanRMS(graph)[0]) 
    slot_data_spe_L_err[slot].append(GetMeanRMS(graph2)[1]) 
    slot_data_spe_R_err[slot].append(GetMeanRMS(graph)[1]) 

    graph = rootfile.Get('g_lyso_L_pc_per_kev_vs_bar')
    c = ROOT.TCanvas('c', 'c', 800, 800)
    graph.Draw()
    #c.SaveAs(plotdir+str(params[module][0])+"_graph.png")
    for point in range(graph.GetN()):
        h_spe_L_ch.Fill(graph.GetPointY(point))
        g_src_vs_LO_ch_L.SetPointX(num*16+point, graph.GetPointY(point))
        #g_spe_vs_asymm_ch_L.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_src_ch_L.SetPointY(num*16+point, graph.GetPointY(point))


    graph = rootfile.Get('g_lyso_R_pc_per_kev_vs_bar'); graph2 = rootfile.Get('g_lyso_L_pc_per_kev_vs_bar')
    for point in range(graph.GetN()):
        h_spe_R_ch.Fill(graph.GetPointY(point))
        g_src_vs_LO_ch_R.SetPointX(num*16+point, graph.GetPointY(point))
        #g_spe_vs_asymm_ch_R.SetPointX(num*16+point, graph.GetPointY(point))
        g_spe_vs_src_ch_R.SetPointY(num*16+point, graph.GetPointY(point))
        g_src_asymm_vs_LO_asymm.SetPointX(num*16+point, 2*(graph2.GetPointY(point)-graph.GetPointY(point))/(graph2.GetPointY(point)+graph.GetPointY(point)))

    
    slot_data_src_L[slot].append(GetMeanRMS(graph2)[0]) 
    slot_data_src_R[slot].append(GetMeanRMS(graph)[0]) 
    slot_data_src_L_err[slot].append(GetMeanRMS(graph2)[1]) 
    slot_data_src_R_err[slot].append(GetMeanRMS(graph)[1]) 
    



    graph = rootfile.Get('g_avg_light_yield_vs_bar')
    h_LO_avg_bar.Fill(GetMeanRMS(graph)[0])
    h_LOrms_bar.Fill(GetMeanRMS(graph)[1]/GetMeanRMS(graph)[0]*100.)
    h_LOmaxvar_bar.Fill(GetMaxVar(graph)/GetMeanRMS(graph)[0]*100.)
    
    slot_data_LO[slot].append(GetMeanRMS(graph)[0])
    slot_data_LO_err[slot].append(GetMeanRMS(graph)[1])

    for point in range(graph.GetN()):
        h_LO_avg_ch.Fill(graph.GetPointY(point))
    
    graph = rootfile.Get('g_light_yield_asymm_vs_bar')
    h_LO_asymm_bar.Fill(GetMeanRMS_abs(graph)[0])
    for point in range(graph.GetN()):
        h_LO_asymm_ch.Fill(graph.GetPointY(point))
        g_spe_vs_asymm_ch_L.SetPointY(num*16+point, graph.GetPointY(point))
        g_spe_vs_asymm_ch_R.SetPointY(num*16+point, graph.GetPointY(point))
        g_LO_vs_asymm_ch_L.SetPointY(num*16+point, graph.GetPointY(point))
        g_LO_vs_asymm_ch_R.SetPointY(num*16+point, graph.GetPointY(point))
        g_spe_asymm_vs_LO_asymm.SetPointY(num*16+point, graph.GetPointY(point))
        g_src_asymm_vs_LO_asymm.SetPointY(num*16+point, graph.GetPointY(point))
            
    graph = rootfile.Get('g_L_light_yield_vs_bar')
    h_LO_L_bar.Fill(GetMeanRMS(graph)[0])
    for point in range(graph.GetN()):
        h_LO_L_ch.Fill(graph.GetPointY(point))
        g_spe_vs_LO_ch_L.SetPointY(num*16+point, graph.GetPointY(point))
        g_LO_vs_asymm_ch_L.SetPointX(num*16+point, graph.GetPointX(point))
        g_src_vs_LO_ch_L.SetPointY(num*16+point, graph.GetPointY(point))

    graph = rootfile.Get('g_R_light_yield_vs_bar')
    h_LO_R_bar.Fill(GetMeanRMS(graph)[0])
    for point in range(graph.GetN()):
        h_LO_R_ch.Fill(graph.GetPointY(point))
        g_spe_vs_LO_ch_R.SetPointY(num*16+point, graph.GetPointY(point))
        g_LO_vs_asymm_ch_R.SetPointX(num*16+point, graph.GetPointX(point))
        g_src_vs_LO_ch_R.SetPointY(num*16+point, graph.GetPointY(point))

    graph = rootfile.Get('g_light_yield_vs_ch')
    h_LOrms_ch.Fill(GetMeanRMS(graph)[1]/GetMeanRMS(graph)[0]*100.)
    h_LOmaxvar_ch.Fill(GetMaxVar(graph)/GetMeanRMS(graph)[0]*100.)

for i in range(12):
    g_slot_vs_spe_L.SetPointX(i,i)
    g_slot_vs_spe_L.SetPointY(i, np.average(np.array(slot_data_spe_L[i])))
    g_slot_vs_spe_L.SetPointError(i, 0, np.std(np.array(slot_data_spe_L[i])))
    g_slot_vs_spe_R.SetPointX(i,i)
    g_slot_vs_spe_R.SetPointY(i, np.average(np.array(slot_data_spe_R[i])))
    g_slot_vs_spe_R.SetPointError(i, 0, np.std(np.array(slot_data_spe_R[i])))
    
    g_slot_vs_src_L.SetPointX(i,i)
    g_slot_vs_src_L.SetPointY(i, np.average(np.array(slot_data_src_L[i])))
    g_slot_vs_src_L.SetPointError(i, 0, np.std(np.array(slot_data_src_L[i])))
    
    g_slot_vs_src_R.SetPointX(i,i)
    g_slot_vs_src_R.SetPointY(i, np.average(np.array(slot_data_src_R[i])))
    g_slot_vs_src_R.SetPointError(i, 0, np.std(np.array(slot_data_src_R[i])))
    
    g_slot_vs_LO.SetPointX(i,i)
    g_slot_vs_LO.SetPointY(i, np.average(np.array(slot_data_LO[i])))
    g_slot_vs_LO.SetPointError(i, 0, np.std(np.array(slot_data_LO[i])))


# draw histos

c = ROOT.TCanvas('c_spe_LR_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_spe_L_ch.SetTitle(';single p.e. charge [pC];entries')
h_spe_L_ch.SetFillStyle(3001)
h_spe_L_ch.SetFillColor(ROOT.kRed)
h_spe_L_ch.SetLineColor(ROOT.kRed)
h_spe_L_ch.GetYaxis().SetRangeUser(0.,1.1*max(h_spe_L_ch.GetMaximum(),h_spe_R_ch.GetMaximum()))
h_spe_L_ch.Draw()
latex_L = ROOT.TLatex(0.64,0.70,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_spe_L_ch.GetMean(),h_spe_L_ch.GetRMS()/h_spe_L_ch.GetMean()*100.))
latex_L.SetNDC()
latex_L.SetTextSize(0.05)
latex_L.SetTextColor(ROOT.kRed)
latex_L.Draw('same')
h_spe_R_ch.SetFillStyle(3001)
h_spe_R_ch.SetFillColor(ROOT.kBlue)
h_spe_R_ch.SetLineColor(ROOT.kBlue)
h_spe_R_ch.Draw('same')
latex_R = ROOT.TLatex(0.64,0.40,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_spe_R_ch.GetMean(),h_spe_R_ch.GetRMS()/h_spe_R_ch.GetMean()*100.))
latex_R.SetNDC()
latex_R.SetTextSize(0.05)
latex_R.SetTextColor(ROOT.kBlue)
latex_R.Draw('same')
c.Print('%s/h_spe_LR_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_avg_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_avg_bar.SetTitle(';avg. bar light output [pe/MeV];entries')
h_LO_avg_bar.SetFillStyle(3001)
h_LO_avg_bar.SetFillColor(ROOT.kBlack)
h_LO_avg_bar.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_avg_bar.GetMean(),h_LO_avg_bar.GetRMS()/h_LO_avg_bar.GetMean()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
line_low = ROOT.TLine(0.85*3200.,0.,0.85*3200,1.05*h_LO_avg_bar.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
c.Print('%s/h_LO_avg_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_avg_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_avg_ch.SetTitle(';bar light output [pe/MeV];entries')
h_LO_avg_ch.SetFillStyle(3001)
h_LO_avg_ch.SetFillColor(ROOT.kBlack)
h_LO_avg_ch.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_avg_ch.GetMean(),h_LO_avg_ch.GetRMS()/h_LO_avg_ch.GetMean()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same')
line_low = ROOT.TLine(0.85*3200.,0.,0.85*3200,1.05*h_LO_avg_ch.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
c.Print('%s/h_LO_avg_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_LR_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_L_bar.SetTitle(';avg. channel light output [pe/MeV];entries')
h_LO_L_bar.SetFillStyle(3001)
h_LO_L_bar.SetFillColor(ROOT.kRed)
h_LO_L_bar.SetLineColor(ROOT.kRed)
h_LO_L_bar.GetYaxis().SetRangeUser(0.,1.1*max(h_LO_L_bar.GetMaximum(),h_LO_R_bar.GetMaximum()))
h_LO_L_bar.Draw()
latex_L = ROOT.TLatex(0.64,0.70,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_L_bar.GetMean(),h_LO_L_bar.GetRMS()/h_LO_L_bar.GetMean()*100.))
latex_L.SetNDC()
latex_L.SetTextSize(0.05)
latex_L.SetTextColor(ROOT.kRed)
latex_L.Draw('same')
h_LO_R_bar.SetFillStyle(3001)
h_LO_R_bar.SetFillColor(ROOT.kBlue)
h_LO_R_bar.SetLineColor(ROOT.kBlue)
h_LO_R_bar.Draw('same')
latex_R = ROOT.TLatex(0.64,0.40,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_R_bar.GetMean(),h_LO_R_bar.GetRMS()/h_LO_R_bar.GetMean()*100.))
latex_R.SetNDC()
latex_R.SetTextSize(0.05)
latex_R.SetTextColor(ROOT.kBlue)
latex_R.Draw('same')
line_low = ROOT.TLine(0.85*3200.,0.,0.85*3200.,1.*max(h_LO_L_bar.GetMaximum(),h_LO_R_bar.GetMaximum()))
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
c.Print('%s/h_LO_LR_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_LR_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_L_ch.SetTitle(';channel light output [pe/MeV];entries')
h_LO_L_ch.SetFillStyle(3001)
h_LO_L_ch.SetFillColor(ROOT.kRed)
h_LO_L_ch.SetLineColor(ROOT.kRed)
h_LO_L_ch.GetYaxis().SetRangeUser(0.,1.1*max(h_LO_L_ch.GetMaximum(),h_LO_R_ch.GetMaximum()))
h_LO_L_ch.Draw()
latex_L = ROOT.TLatex(0.64,0.70,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_L_ch.GetMean(),h_LO_L_ch.GetRMS()/h_LO_L_ch.GetMean()*100.))
latex_L.SetNDC()
latex_L.SetTextSize(0.05)
latex_L.SetTextColor(ROOT.kRed)
latex_L.Draw('same')
h_LO_R_ch.SetFillStyle(3001)
h_LO_R_ch.SetFillColor(ROOT.kBlue)
h_LO_R_ch.SetLineColor(ROOT.kBlue)
h_LO_R_ch.Draw('same')
latex_R = ROOT.TLatex(0.64,0.40,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_R_ch.GetMean(),h_LO_R_ch.GetRMS()/h_LO_R_ch.GetMean()*100.))
latex_R.SetNDC()
latex_R.SetTextSize(0.05)
latex_R.SetTextColor(ROOT.kBlue)
latex_R.Draw('same')
line_low = ROOT.TLine(0.85*3200.,0.,0.85*3200.,1.*max(h_LO_L_ch.GetMaximum(),h_LO_R_ch.GetMaximum()))
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
c.Print('%s/h_LO_LR_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_asymm_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_asymm_bar.SetTitle(';avg. L.O. asymmetry [ 2*(L-R)/(L+R) ];entries')
h_LO_asymm_bar.SetFillStyle(3001)
h_LO_asymm_bar.SetFillColor(ROOT.kBlack)
h_LO_asymm_bar.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_asymm_bar.GetMean(),h_LO_asymm_bar.GetRMS()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
# line_low = ROOT.TLine(-0.1,0.,-0.1,1.05*h_LO_asymm_bar.GetMaximum())
# line_low.SetLineColor(ROOT.kGreen+1)
# line_low.SetLineWidth(4)
# line_low.SetLineStyle(2)
# line_low.Draw('same')
line_high = ROOT.TLine(0.1,0.,0.1,1.05*h_LO_asymm_bar.GetMaximum())
line_high.SetLineColor(ROOT.kGreen+1)
line_high.SetLineWidth(4)
line_high.SetLineStyle(2)
line_high.Draw('same')
c.Print('%s/h_LO_asymm_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_asymm_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LO_asymm_ch.SetTitle(';L.O. asymmetry [ 2*(L-R)/(L+R) ];entries')
h_LO_asymm_ch.SetFillStyle(3001)
h_LO_asymm_ch.SetFillColor(ROOT.kBlack)
h_LO_asymm_ch.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_asymm_ch.GetMean(),h_LO_asymm_ch.GetRMS()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
line_low = ROOT.TLine(-0.2,0.,-0.2,1.05*h_LO_asymm_ch.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
line_high = ROOT.TLine(0.2,0.,0.2,1.05*h_LO_asymm_ch.GetMaximum())
line_high.SetLineColor(ROOT.kGreen+1)
line_high.SetLineWidth(4)
line_high.SetLineStyle(2)
line_high.Draw('same')
c.Print('%s/h_LO_asymm_ch.png'%plotDir)




c = ROOT.TCanvas('c_LOrms_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LOrms_bar.SetTitle(';bar RMS [%];entries')
h_LOrms_bar.SetFillStyle(3001)
h_LOrms_bar.SetFillColor(ROOT.kBlack)
h_LOrms_bar.Draw()
line = ROOT.TLine(5.,0.,5.,1.05*h_LOrms_bar.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
c.Print('%s/h_LOrms_bar.png'%plotDir)

c = ROOT.TCanvas('c_LOrms_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LOrms_ch.SetTitle(';channel RMS [%];entries')
h_LOrms_ch.SetFillStyle(3001)
h_LOrms_ch.SetFillColor(ROOT.kBlack)
h_LOrms_ch.Draw()
line = ROOT.TLine(7.,0.,7.,1.05*h_LOrms_ch.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
c.Print('%s/h_LOrms_ch.png'%plotDir)



c = ROOT.TCanvas('c_LOmaxvar_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LOmaxvar_bar.SetTitle(';bar max. var. [%];entries')
h_LOmaxvar_bar.SetFillStyle(3001)
h_LOmaxvar_bar.SetFillColor(ROOT.kBlack)
h_LOmaxvar_bar.Draw()
line = ROOT.TLine(30.,0.,30.,1.05*h_LOmaxvar_bar.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
c.Print('%s/h_LOmaxvar_bar.png'%plotDir)

c = ROOT.TCanvas('c_LOmaxvar_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_LOmaxvar_ch.SetTitle(';channel max. var. [%];entries')
h_LOmaxvar_ch.SetFillStyle(3001)
h_LOmaxvar_ch.SetFillColor(ROOT.kBlack)
h_LOmaxvar_ch.Draw()
line = ROOT.TLine(40.,0.,40.,1.05*h_LOmaxvar_ch.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
c.Print('%s/h_LOmaxvar_ch.png'%plotDir)

c = ROOT.TCanvas('c_SPE_vs_LO_ch','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_spe_vs_LO_ch_L.SetTitle("'';SPE Charge (pC);Light Output (PE/Mev)")
#g_spe_vs_LO_ch_L.GetXaxis().SetLabelSize(0.35)
#g_spe_vs_LO_ch_L.GetYaxis().SetLabelSize(0.35)
g_spe_vs_LO_ch_L.SetFillStyle(3001)
g_spe_vs_LO_ch_L.SetMarkerColor(ROOT.kRed)
g_spe_vs_LO_ch_L.Draw("ap")
g_spe_vs_LO_ch_L.GetXaxis().SetRangeUser(2.8, 4.2)
g_spe_vs_LO_ch_L.GetYaxis().SetRangeUser(2000, 4000)
#g_spe_vs_LO_ch_R.SetTitle("Light Output (PE/MEV);SPE Charge (pC);''")
g_spe_vs_LO_ch_R.SetFillStyle(3001)
g_spe_vs_LO_ch_R.SetMarkerColor(ROOT.kBlue)
g_spe_vs_LO_ch_R.Draw("p same")
c.SetLeftMargin(0.2)
c.Print('%s/g_spe_vs_LO_ch.png'%plotDir)



c = ROOT.TCanvas('c_SPE_vs_src_ch','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_spe_vs_src_ch_L.SetTitle("'';SPE Charge (pC);SRC Charge (PE/keV)")
g_spe_vs_src_ch_L.SetFillStyle(3001)
g_spe_vs_src_ch_L.SetMarkerColor(ROOT.kRed)
g_spe_vs_src_ch_L.GetXaxis().SetRangeUser(2.8, 4.2)
g_spe_vs_src_ch_L.GetYaxis().SetRangeUser(1.4, 2.4)
g_spe_vs_src_ch_L.Draw("ap")
#g_spe_vs_src_ch_R.SetTitle("Light Output (PE/MEV);SPE Charge (pC);''")
g_spe_vs_src_ch_R.SetFillStyle(3001)
g_spe_vs_src_ch_R.SetMarkerColor(ROOT.kBlue)
g_spe_vs_src_ch_R.Draw("p same")
c.Print('%s/g_spe_vs_src_ch.png'%plotDir)


c = ROOT.TCanvas('c_src_vs_LO_ch','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_src_vs_LO_ch_L.SetTitle("'';SRC Charge (PE/keV);Light Output (PE/MEV)")
g_src_vs_LO_ch_L.SetFillStyle(3001)
g_src_vs_LO_ch_L.SetMarkerColor(ROOT.kRed)
g_src_vs_LO_ch_L.GetXaxis().SetRangeUser(1.4, 2.4)
g_src_vs_LO_ch_L.GetYaxis().SetRangeUser(2000, 4000)
g_src_vs_LO_ch_L.Draw("ap")
#g_src_vs_LO_ch_R.SetTitle("Light Output (PE/MEV);SPE Charge (pC);''")
g_src_vs_LO_ch_R.SetFillStyle(3001)
g_src_vs_LO_ch_R.SetMarkerColor(ROOT.kBlue)
g_src_vs_LO_ch_R.Draw("p same")

c.Print('%s/g_src_vs_LO_ch.png'%plotDir)


c = ROOT.TCanvas('c_spe_asymm_vs_LO_asymm','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_spe_asymm_vs_LO_asymm.SetTitle("'';SPE Asymm ;LO Asymm")
g_spe_asymm_vs_LO_asymm.SetFillStyle(3001)
g_spe_asymm_vs_LO_asymm.SetMarkerColor(ROOT.kBlack)
g_spe_asymm_vs_LO_asymm.Draw("ap")

c.Print('%s/g_spe_asymm_vs_LO_asymm.png'%plotDir)


c = ROOT.TCanvas('c_src_asymm_vs_LO_asymm','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_src_asymm_vs_LO_asymm.SetTitle("'';SRC Asymm ;Asymm")
g_src_asymm_vs_LO_asymm.SetFillStyle(3001)
g_src_asymm_vs_LO_asymm.SetMarkerColor(ROOT.kBlack)
g_src_asymm_vs_LO_asymm.Draw("ap")

c.Print('%s/g_src_asymm_vs_LO_asymm.png'%plotDir)




c = ROOT.TCanvas('c_slot_vs_spe','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_slot_vs_spe_L.SetTitle("'';Slot;SPE Charge (PC)")
g_slot_vs_spe_L.SetFillStyle(3001)
g_slot_vs_spe_L.SetMarkerColor(ROOT.kRed)
g_slot_vs_spe_L.SetLineColor(ROOT.kRed)
g_slot_vs_spe_L.Draw("ap")
#g_src_vs_LO_ch_R.SetTitle("Light Output (PE/MEV);SPE Charge (pC);''")
g_slot_vs_spe_R.SetFillStyle(3001)
g_slot_vs_spe_R.SetMarkerColor(ROOT.kBlue)
g_slot_vs_spe_R.SetLineColor(ROOT.kBlue)
#g_slot_vs_spe_R.GetYaxis().SetMaximum(4)
g_slot_vs_spe_R.Draw("p same")
g_slot_vs_spe_R.GetYaxis().SetRangeUser(3.3,3.8)
g_slot_vs_spe_L.GetYaxis().SetRangeUser(3.3,3.8)


c.Print('%s/g_slot_vs_spe.png'%plotDir)





c = ROOT.TCanvas('c_Slot_vs_LO','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_slot_vs_LO.SetTitle("'';Slot ;Light Output")
g_slot_vs_LO.SetFillStyle(3001)
g_slot_vs_LO.SetMarkerColor(ROOT.kBlack)
g_slot_vs_LO.SetLineColor(ROOT.kBlack)
g_slot_vs_LO.Draw("ap")
g_slot_vs_LO.GetYaxis().SetRangeUser(3000, 3500)



c.Print('%s/g_slot_vs_LO.png'%plotDir)



c = ROOT.TCanvas('c_slot_vs_src','',800,800)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_slot_vs_src_L.SetTitle("'';Slot;SRC Charge (pC/keV)")
g_slot_vs_src_L.SetFillStyle(3001)
g_slot_vs_src_L.SetMarkerColor(ROOT.kRed)
g_slot_vs_src_L.SetLineColor(ROOT.kRed)
g_slot_vs_src_L.Draw("ap")
#g_src_vs_LO_ch_R.SetTitle("Light Output (PE/MEV);SPE Charge (pC);''")
g_slot_vs_src_R.SetFillStyle(3001)
g_slot_vs_src_R.SetMarkerColor(ROOT.kBlue)
g_slot_vs_src_R.SetLineColor(ROOT.kBlue)
#g_slot_vs_src_R.GetYaxis().SetMaximum(4)
g_slot_vs_src_R.Draw("p same")
#g_slot_vs_src_R.GetYaxis().SetRangeUser(3.3,3.8)
#g_slot_vs_src_L.GetYaxis().SetRangeUser(3.3,3.8)


c.Print('%s/g_slot_vs_src.png'%plotDir)

print(f"{len(modules_tested)} modules tested")
