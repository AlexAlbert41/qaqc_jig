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

#from typing import NamedTuple


#data_path = '/data1/SMQAQC/PRODUCTION/'
#selections = []
#plotDir = '/data1/html/data1/SMQAQC/PRODUCTION/summaryPlots_SMID_42to152/'
#plotDir = '/data1/html/data1/SMQAQC/PRODUCTION/summaryPlots_SMID_1to152/'

data_path = '/data/QAQC_SM/qaqc-gui_output/Summary_Plots_First37_NoCal/SM_QAQC_Production/'
selections = []
plotDir = '/data/QAQC_SM/qaqc-gui_output/Summary_Plots_NewDirectory_NoCal_AddedPlots_Ranges/'


#runs = ["63-66","85-88","90-93","95-101","104-104"] # 42-152
#runs = ["50-56","63-66","85-88","90-93","95-101","104-104"] # 1-152
runs = ["10-200"]
#modules_acc = ["42-152"] # 42-152
modules_acc = ["1-152"] # 1-152


MIN_SPE_ch = -1000 # modified by Alex
MAX_SPE_ch = 1000
MIN_LO_bar = 0.90 * 3150.
MIN_LO_ch = 0.85 * 3150.
MAX_LO_ASYMM_bar = 0.06
MIN_LO_ASYMM_ch = -0.15
MAX_LO_ASYMM_ch = 0.15


#MIN_SPE_ch = -1000 # modified by Alex
#MAX_SPE_ch = 1000
#MIN_LO_bar = 0
#MIN_LO_ch = 0
#MAX_LO_ASYMM_bar = 1
#MIN_LO_ASYMM_ch = -0.15
#MAX_LO_ASYMM_ch = 0.15



nCatA = 0
nCatB = 0

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.25,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)

def expand_range(rng):
    numbers = []
    for item in rng:
        if '-' in item:
            start, end = map(int, item.split('-'))   
            numbers.extend(range(start, end + 1))
        else:
            numbers.append(int(item))
    return sorted(numbers)

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

# Create a list of all runs
list_runs = expand_range(runs)
list_modules = expand_range(modules_acc)
# string list
list_runs = list(map(str, list_runs))
list_modules = list(map(str, list_modules))
# add barcode
prefix = "32110020"
list_modules = ["{}{:06d}".format(prefix, int(mod)) for mod in list_modules]
print(list_modules)
# retrieving root files 
inputFiles = glob.glob(data_path+'/run*/*_analysis.root')
for inputFile in inputFiles:
    tokens = inputFile.split('/')
    run = ''
    for token in tokens:
        if 'module' in token:
            module = token[7:21] # SM ID
        if 'run' in token:
            run = int(token[3:]) # run number
    '''
    # saving in the dict only selected runs and modules
    if str(run) in list_runs and str(module) in list_modules:
        if run == 56 and '32110020000030' in module:
            continue
        print("run ", run, "  module ", module)
        #modules.append(module)
        params[module] = [inputFile,run,'GOOD']
    '''
    params[module] = [inputFile,run,'GOOD']

    modules.append(module)
print(modules)
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
h_LO_asymm_ch = ROOT.TH1F('h_LO_asymm_ch','',100,-0.3,0.3)

h_LOrms_bar = ROOT.TH1F('h_LOrms_bar','',60,0.,30.)
h_LOrms_ch = ROOT.TH1F('h_LOrms_ch','',60,0.,30.)

h_LOmaxvar_bar = ROOT.TH1F('h_LOmaxvar_bar','',50,0.,100.)
h_LOmaxvar_ch = ROOT.TH1F('h_LOmaxvar_ch','',50,0.,100.)

modules_tested = ["32110020000016"]

# selecting the modules to be included in the summary: accept 1 if included, 0 otherwise
for module in modules:
    print(module)
    if module in modules_tested:
        continue
    modules_tested.append(module)
    param = params[module]
    # accept = 1
    # if param[1] < 50: # measurements re-done after changing module boards
    #     accept = 0 
    # elif param[1] > 56 and param[1] < 63: # uniformity tests
    #     accept = 0 
    # elif param[1]>67:
    #     accept = 0
    
    # #for selection in selections:
    # #    tempAccept = 0
    # #    for param in params:
    # #        if selection in param:
    # #            tempAccept = 1
    # #    accept *= tempAccept
    # if accept == 0:
    #     continue
    
    print("analyzing ", module,"    --> ",param)
    
    
    rootfile = ROOT.TFile(params[module][0],'READ')
    
    # acceptance
    isCatA = True
    isCatB = False
    
    # filling histos
    graph = rootfile.Get('g_spe_L_vs_bar')
    for point in range(graph.GetN()):
        h_spe_L_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_SPE_ch or graph.GetPointY(point) > MAX_SPE_ch:
            isCatB = True
            isCatA = False
    
    graph = rootfile.Get('g_spe_R_vs_bar')
    for point in range(graph.GetN()):
        h_spe_R_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_SPE_ch or graph.GetPointY(point) > MAX_SPE_ch:
            isCatB = True
            isCatA = False
    
    graph = rootfile.Get('g_avg_light_yield_vs_bar')
    h_LO_avg_bar.Fill(GetMeanRMS(graph)[0])
    h_LOrms_bar.Fill(GetMeanRMS(graph)[1]/GetMeanRMS(graph)[0]*100.)
    h_LOmaxvar_bar.Fill(GetMaxVar(graph)/GetMeanRMS(graph)[0]*100.)
    h_LO_avg_ch.Fill(graph.GetPointY(point))
    if GetMeanRMS(graph)[0] < MIN_LO_bar:
        isCatB = True
        isCatA = False
    for point in range(graph.GetN()):
        h_LO_avg_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_LO_bar:
            isCatB = True
            isCatA = False
        
    graph = rootfile.Get('g_light_yield_asymm_vs_bar')
    h_LO_asymm_bar.Fill(GetMeanRMS_abs(graph)[0])
    if GetMeanRMS_abs(graph)[0] > MAX_LO_ASYMM_bar:
        isCatB = True
        isCatA = False
    
    for point in range(graph.GetN()):
        h_LO_asymm_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_LO_ASYMM_ch or graph.GetPointY(point) > MAX_LO_ASYMM_ch:
            isCatB = True
            isCatA = False
            
    graph = rootfile.Get('g_L_light_yield_vs_bar')
    h_LO_L_bar.Fill(GetMeanRMS(graph)[0])
    if GetMeanRMS(graph)[0] < MIN_LO_ch:
        isCatB = True
        isCatA = False
    for point in range(graph.GetN()):
        h_LO_L_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_LO_ch:
            isCatB = True
            isCatA = False
    
    graph = rootfile.Get('g_R_light_yield_vs_bar')
    h_LO_R_bar.Fill(GetMeanRMS(graph)[0])
    if GetMeanRMS(graph)[0] < MIN_LO_ch:
        isCatB = True
        isCatA = False
    for point in range(graph.GetN()):
        h_LO_R_ch.Fill(graph.GetPointY(point))
        if graph.GetPointY(point) < MIN_LO_ch:
            isCatB = True
            isCatA = False
    
    graph = rootfile.Get('g_light_yield_vs_ch')
    h_LOrms_ch.Fill(GetMeanRMS(graph)[1]/GetMeanRMS(graph)[0]*100.)
    h_LOmaxvar_ch.Fill(GetMaxVar(graph)/GetMeanRMS(graph)[0]*100.)
    

    if isCatB == isCatA:
        print('Error: a module cannot be both catA and catB')
    if isCatA:
        nCatA += 1
    if isCatB:
        nCatB += 1
        print('>>> '+module+' is cat. B!!!')

# draw histos

latex_cat = ROOT.TLatex(0.18,0.85,'#splitline{cat. A: %d (%.1f%%)}{cat. B: %d (%.1f%%)}'%(nCatA,100.*nCatA/(nCatA+nCatB),nCatB,100.*nCatB/(nCatA+nCatB)))
latex_cat.SetNDC()
latex_cat.SetTextSize(0.05)

c = ROOT.TCanvas('c_spe_LR_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_spe_L_ch.SetTitle(';single p.e. charge [pC];entries')
h_spe_L_ch.SetFillStyle(3001)
h_spe_L_ch.SetFillColor(ROOT.kRed)
h_spe_L_ch.SetLineColor(ROOT.kRed)
h_spe_L_ch.GetYaxis().SetRangeUser(0.5,1.1*max(h_spe_L_ch.GetMaximum(),h_spe_R_ch.GetMaximum()))
h_spe_L_ch.GetXaxis().SetRangeUser(2, 4.2)
h_spe_L_ch.Draw()
latex_L = ROOT.TLatex(0.64,0.70,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_spe_L_ch.GetMean(),h_spe_L_ch.GetRMS()/h_spe_L_ch.GetMean()*100.))
latex_L.SetNDC()
latex_L.SetTextSize(0.05)
latex_L.SetTextColor(ROOT.kRed)
latex_L.Draw('same')
h_spe_R_ch.SetFillStyle(3001)
h_spe_R_ch.SetFillColor(ROOT.kBlue)
h_spe_R_ch.GetXaxis().SetRangeUser(2.8, 4.2)
h_spe_R_ch.SetLineColor(ROOT.kBlue)
h_spe_R_ch.Draw('same')
latex_R = ROOT.TLatex(0.64,0.40,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_spe_R_ch.GetMean(),h_spe_R_ch.GetRMS()/h_spe_R_ch.GetMean()*100.))
latex_R.SetNDC()
latex_R.SetTextSize(0.05)
latex_R.SetTextColor(ROOT.kBlue)
latex_R.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_spe_LR_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_avg_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_avg_bar.SetTitle(';avg. bar light output [pe/MeV];entries')
h_LO_avg_bar.SetFillStyle(3001)
h_LO_avg_bar.SetFillColor(ROOT.kBlack)
h_LO_avg_bar.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_avg_bar.GetMean(),h_LO_avg_bar.GetRMS()/h_LO_avg_bar.GetMean()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
line_low = ROOT.TLine(MIN_LO_bar,0.,MIN_LO_bar,1.05*h_LO_avg_bar.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_avg_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_avg_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_avg_ch.SetTitle(';bar light output [pe/MeV];entries')
h_LO_avg_ch.SetFillStyle(3001)
h_LO_avg_ch.SetFillColor(ROOT.kBlack)
h_LO_avg_ch.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_avg_ch.GetMean(),h_LO_avg_ch.GetRMS()/h_LO_avg_ch.GetMean()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same')
line_low = ROOT.TLine(MIN_LO_bar,0.,MIN_LO_bar,1.05*h_LO_avg_ch.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_avg_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_LR_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_L_bar.SetTitle(';avg. channel light output [pe/MeV];entries')
h_LO_L_bar.SetFillStyle(3001)
h_LO_L_bar.SetFillColor(ROOT.kRed)
h_LO_L_bar.SetLineColor(ROOT.kRed)
h_LO_L_bar.GetYaxis().SetRangeUser(0.5,1.1*max(h_LO_L_bar.GetMaximum(),h_LO_R_bar.GetMaximum()))
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
line_low = ROOT.TLine(MIN_LO_ch,0.,MIN_LO_ch,1.*max(h_LO_L_bar.GetMaximum(),h_LO_R_bar.GetMaximum()))
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_LR_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_LR_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_L_ch.SetTitle(';channel light output [pe/MeV];entries')
h_LO_L_ch.SetFillStyle(3001)
h_LO_L_ch.SetFillColor(ROOT.kRed)
h_LO_L_ch.SetLineColor(ROOT.kRed)
h_LO_L_ch.GetYaxis().SetRangeUser(0.5,1.1*max(h_LO_L_ch.GetMaximum(),h_LO_R_ch.GetMaximum()))
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
line_low = ROOT.TLine(MIN_LO_ch,0.,MIN_LO_ch,1.*max(h_LO_L_ch.GetMaximum(),h_LO_R_ch.GetMaximum()))
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_LR_ch.png'%plotDir)




c = ROOT.TCanvas('c_LO_asymm_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_asymm_bar.SetTitle(';avg. L.O. asymmetry [ 2*(L-R)/(L+R) ];entries')
h_LO_asymm_bar.SetFillStyle(3001)
h_LO_asymm_bar.SetFillColor(ROOT.kBlack)
h_LO_asymm_bar.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_asymm_bar.GetMean(),h_LO_asymm_bar.GetRMS()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
line_high = ROOT.TLine(MAX_LO_ASYMM_bar,0.,MAX_LO_ASYMM_bar,1.05*h_LO_asymm_bar.GetMaximum())
line_high.SetLineColor(ROOT.kGreen+1)
line_high.SetLineWidth(4)
line_high.SetLineStyle(2)
line_high.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_asymm_bar.png'%plotDir)

c = ROOT.TCanvas('c_LO_asymm_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LO_asymm_ch.SetTitle(';L.O. asymmetry [ 2*(L-R)/(L+R) ];entries')
h_LO_asymm_ch.SetFillStyle(3001)
h_LO_asymm_ch.SetFillColor(ROOT.kBlack)
h_LO_asymm_ch.Draw()
latex = ROOT.TLatex(0.64,0.60,'#splitline{mean: %.2e}{RMS: %.1f %%}'%(h_LO_asymm_ch.GetMean(),h_LO_asymm_ch.GetRMS()*100.))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.Draw('same') 
line_low = ROOT.TLine(MIN_LO_ASYMM_ch,0.,MIN_LO_ASYMM_ch,1.05*h_LO_asymm_ch.GetMaximum())
line_low.SetLineColor(ROOT.kGreen+1)
line_low.SetLineWidth(4)
line_low.SetLineStyle(2)
line_low.Draw('same')
line_high = ROOT.TLine(MAX_LO_ASYMM_ch,0.,MAX_LO_ASYMM_ch,1.05*h_LO_asymm_ch.GetMaximum())
line_high.SetLineColor(ROOT.kGreen+1)
line_high.SetLineWidth(4)
line_high.SetLineStyle(2)
line_high.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LO_asymm_ch.png'%plotDir)




c = ROOT.TCanvas('c_LOrms_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LOrms_bar.SetTitle(';bar RMS [%];entries')
h_LOrms_bar.SetFillStyle(3001)
h_LOrms_bar.SetFillColor(ROOT.kBlack)
h_LOrms_bar.Draw()
line = ROOT.TLine(5.,0.,5.,1.05*h_LOrms_bar.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LOrms_bar.png'%plotDir)

c = ROOT.TCanvas('c_LOrms_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LOrms_ch.SetTitle(';channel RMS [%];entries')
h_LOrms_ch.SetFillStyle(3001)
h_LOrms_ch.SetFillColor(ROOT.kBlack)
h_LOrms_ch.Draw()
line = ROOT.TLine(7.,0.,7.,1.05*h_LOrms_ch.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LOrms_ch.png'%plotDir)



c = ROOT.TCanvas('c_LOmaxvar_bar','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LOmaxvar_bar.SetTitle(';bar max. var. [%];entries')
h_LOmaxvar_bar.SetFillStyle(3001)
h_LOmaxvar_bar.SetFillColor(ROOT.kBlack)
h_LOmaxvar_bar.Draw()
line = ROOT.TLine(30.,0.,30.,1.05*h_LOmaxvar_bar.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LOmaxvar_bar.png'%plotDir)

c = ROOT.TCanvas('c_LOmaxvar_ch','',800,700)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetLogy()
h_LOmaxvar_ch.SetTitle(';channel max. var. [%];entries')
h_LOmaxvar_ch.SetFillStyle(3001)
h_LOmaxvar_ch.SetFillColor(ROOT.kBlack)
h_LOmaxvar_ch.Draw()
line = ROOT.TLine(40.,0.,40.,1.05*h_LOmaxvar_ch.GetMaximum())
line.SetLineColor(ROOT.kGreen+1)
line.SetLineWidth(4)
line.SetLineStyle(2)
line.Draw('same')
latex_cat.Draw('same')
c.Print('%s/h_LOmaxvar_ch.png'%plotDir)
