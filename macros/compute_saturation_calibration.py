#! /usr/bin/env python3                                                                                                                                                                                                                                                         
import os
import shutil
import glob
import math
import array
import sys
import time
import json
import numpy as np
from pathlib import Path

import ROOT
import tdrstyle


sodium_charge = 511
cesium_charge = 662

input_directory = Path("/run/media/cptlab/qaqc_data/results/QAQC_SM/qaqc-gui_output")
plotDir = "saturation_correction_plots/"
ROOT.gROOT.SetBatch(True)

if not os.path.isdir(plotDir):os.mkdir(plotDir)



#select uncalibrated root files path corresponding to the same module
root_file_paths = list(input_directory.glob(f"**/module*analysis.root"))
modules = []
modules_path_Cs, modules_path_Na_2OV_oldDig, modules_path_Na_2OV_newDig, modules_path_Na_3OV_oldDig, modules_path_Na_3OV_newDig = {},{},{},{},{}
dict_list = [modules_path_Cs, modules_path_Na_2OV_oldDig, modules_path_Na_2OV_newDig, modules_path_Na_3OV_oldDig, modules_path_Na_3OV_newDig]
for file_path in root_file_paths:
    inputFile = str(file_path)
    module_ID = int(inputFile[inputFile.rfind("module_")+len("module_"):inputFile.find("_analysis.root")])
    #print(inputFile)
    if module_ID not in modules:
        modules.append(module_ID)
    #if module_ID!=32110020008593: ##COMMENT OUT IF TESTING MORE MODULES, AND CHANGE OUTPUT PATH
    #    continue
    #print(module_ID)
    if "Na_2Ov_AdjustDig" in inputFile:
        modules_path_Na_2OV_newDig[module_ID] = inputFile
    if "Na_2Ov_NormalDig" in inputFile:
        modules_path_Na_2OV_oldDig[module_ID] = inputFile
    if "Na_3Ov_AdjustDig" in inputFile:
        modules_path_Na_3OV_newDig[module_ID] = inputFile
    if "Na_3Ov_NormalDig" in inputFile:
        modules_path_Na_3OV_oldDig[module_ID] = inputFile
    if "Cs_3Ov_NormalDig" in inputFile:
        modules_path_Cs[module_ID] = inputFile

print("Modules: ", modules)

Na_Graph = ROOT.TGraphErrors()
Na_Graph_LO = ROOT.TGraphErrors()
Na_Graph_SPE = ROOT.TGraphErrors()

#get ratios of source charge for 3 OV sodium first
counter=0
for module, satfile in modules_path_Na_3OV_oldDig.items():
    unsatfile = modules_path_Na_3OV_newDig[module]
    rootfile_sat = ROOT.TFile(satfile,'READ')
    rootfile_unsat = ROOT.TFile(unsatfile,'READ')
    satGraph = rootfile_sat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    unsatGraph = rootfile_unsat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    satGraph_LO = rootfile_sat.Get("g_light_yield_vs_ch")
    unsatGraph_LO = rootfile_unsat.Get("g_light_yield_vs_ch")
    satGraph_SPE = rootfile_sat.Get("g_spe_vs_ch")
    unsatGraph_SPE = rootfile_unsat.Get("g_spe_vs_ch")

    for point in range(satGraph.GetN()):
        Na_Graph.SetPointX(point+32*counter, satGraph.GetPointY(point))
        Na_Graph.SetPointY(point+32*counter, unsatGraph.GetPointY(point)/satGraph.GetPointY(point))
        Na_Graph.SetPointError(point+32*counter, satGraph.GetErrorY(point), unsatGraph.GetPointY(point)/satGraph.GetPointY(point)*np.sqrt((satGraph.GetErrorY(point)/satGraph.GetPointY(point))**2+(unsatGraph.GetErrorY(point)/unsatGraph.GetPointY(point))**2))
        Na_Graph_LO.SetPointX(point+32*counter, satGraph_LO.GetPointY(point))
        Na_Graph_LO.SetPointY(point+32*counter, unsatGraph_LO.GetPointY(point)/satGraph_LO.GetPointY(point))
        Na_Graph_LO.SetPointError(point+32*counter, satGraph_LO.GetErrorY(point), unsatGraph_LO.GetPointY(point)/satGraph_LO.GetPointY(point)*np.sqrt((satGraph_LO.GetErrorY(point)/satGraph_LO.GetPointY(point))**2+(unsatGraph_LO.GetErrorY(point)/unsatGraph_LO.GetPointY(point))**2))
        Na_Graph_SPE.SetPointX(point+32*counter, satGraph_SPE.GetPointY(point))
        Na_Graph_SPE.SetPointY(point+32*counter, unsatGraph_SPE.GetPointY(point)/satGraph_SPE.GetPointY(point))
        Na_Graph_SPE.SetPointError(point+32*counter, satGraph_SPE.GetErrorY(point), unsatGraph_SPE.GetPointY(point)/satGraph_SPE.GetPointY(point)*np.sqrt((satGraph_SPE.GetErrorY(point)/satGraph_SPE.GetPointY(point))**2+(unsatGraph_SPE.GetErrorY(point)/unsatGraph_SPE.GetPointY(point))**2))
    counter+=1

Na_Graph_2OV = ROOT.TGraphErrors()
Na_Graph_LO_2OV = ROOT.TGraphErrors()
Na_Graph_SPE_2OV = ROOT.TGraphErrors()

#get ratios of source charge for 2 OV sodium
counter=0
for module, satfile in modules_path_Na_2OV_oldDig.items():
    unsatfile = modules_path_Na_2OV_newDig[module]
    rootfile_sat = ROOT.TFile(satfile,'READ')
    rootfile_unsat = ROOT.TFile(unsatfile,'READ')
    satGraph = rootfile_sat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    unsatGraph = rootfile_unsat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    satGraph_LO = rootfile_sat.Get("g_light_yield_vs_ch")
    unsatGraph_LO = rootfile_unsat.Get("g_light_yield_vs_ch")
    satGraph_SPE = rootfile_sat.Get("g_spe_vs_ch")
    unsatGraph_SPE = rootfile_unsat.Get("g_spe_vs_ch")
    for point in range(satGraph.GetN()):
        Na_Graph_2OV.SetPointX(point+32*counter, satGraph.GetPointY(point))
        Na_Graph_2OV.SetPointY(point+32*counter, unsatGraph.GetPointY(point)/satGraph.GetPointY(point))
        Na_Graph_2OV.SetPointError(point+32*counter, satGraph.GetErrorY(point), unsatGraph.GetPointY(point)/satGraph.GetPointY(point)*np.sqrt((satGraph.GetErrorY(point)/satGraph.GetPointY(point))**2+(unsatGraph.GetErrorY(point)/unsatGraph.GetPointY(point))**2))
        Na_Graph_LO_2OV.SetPointX(point+32*counter, satGraph_LO.GetPointY(point))
        Na_Graph_LO_2OV.SetPointY(point+32*counter, unsatGraph_LO.GetPointY(point)/satGraph_LO.GetPointY(point))
        Na_Graph_LO_2OV.SetPointError(point+32*counter, satGraph_LO.GetErrorY(point), unsatGraph_LO.GetPointY(point)/satGraph_LO.GetPointY(point)*np.sqrt((satGraph_LO.GetErrorY(point)/satGraph_LO.GetPointY(point))**2+(unsatGraph_LO.GetErrorY(point)/unsatGraph_LO.GetPointY(point))**2))
        Na_Graph_SPE_2OV.SetPointX(point+32*counter, satGraph_SPE.GetPointY(point))
        Na_Graph_SPE_2OV.SetPointY(point+32*counter, unsatGraph_SPE.GetPointY(point)/satGraph_SPE.GetPointY(point))
        Na_Graph_SPE_2OV.SetPointError(point+32*counter, satGraph_SPE.GetErrorY(point), unsatGraph_SPE.GetPointY(point)/satGraph_SPE.GetPointY(point)*np.sqrt((satGraph_SPE.GetErrorY(point)/satGraph_SPE.GetPointY(point))**2+(unsatGraph_SPE.GetErrorY(point)/unsatGraph_SPE.GetPointY(point))**2))
    counter+=1



Cs_Graph_3OV = ROOT.TGraphErrors()
Cs_Graph_LO_3OV = ROOT.TGraphErrors()
Cs_Graph_SPE_3OV = ROOT.TGraphErrors()

#get ratios of source charge for 3 OV cesium with 3 Ov sodium, new digi (to scale data that was already taken)
counter=0
for module, satfile in modules_path_Cs.items():
    unsatfile = modules_path_Na_3OV_newDig[module]
    rootfile_sat = ROOT.TFile(satfile,'READ')
    rootfile_unsat = ROOT.TFile(unsatfile,'READ')
    satGraph = rootfile_sat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    unsatGraph = rootfile_unsat.Get('g_lyso_pc_per_kev_raw_vs_ch')
    satGraph_LO = rootfile_sat.Get("g_light_yield_vs_ch")
    unsatGraph_LO = rootfile_unsat.Get("g_light_yield_vs_ch")
    satGraph_SPE = rootfile_sat.Get("g_spe_vs_ch")
    unsatGraph_SPE = rootfile_unsat.Get("g_spe_vs_ch")
    for point in range(satGraph.GetN()):
        Cs_Graph_3OV.SetPointX(point+32*counter, satGraph.GetPointY(point))
        Cs_Graph_3OV.SetPointY(point+32*counter, unsatGraph.GetPointY(point)/(satGraph.GetPointY(point)))
        Cs_Graph_3OV.SetPointError(point+32*counter, satGraph.GetErrorY(point), unsatGraph.GetPointY(point)/(satGraph.GetPointY(point))*np.sqrt((satGraph.GetErrorY(point)/satGraph.GetPointY(point))**2+(unsatGraph.GetErrorY(point)/unsatGraph.GetPointY(point))**2))
        Cs_Graph_LO_3OV.SetPointX(point+32*counter, satGraph_LO.GetPointY(point))
        Cs_Graph_LO_3OV.SetPointY(point+32*counter, unsatGraph_LO.GetPointY(point)/(satGraph_LO.GetPointY(point)))
        Cs_Graph_LO_3OV.SetPointError(point+32*counter, satGraph_LO.GetErrorY(point), unsatGraph_LO.GetPointY(point)/(satGraph_LO.GetPointY(point))*np.sqrt((satGraph_LO.GetErrorY(point)/satGraph_LO.GetPointY(point))**2+(unsatGraph_LO.GetErrorY(point)/unsatGraph_LO.GetPointY(point))**2))
        Cs_Graph_SPE_3OV.SetPointX(point+32*counter, satGraph_SPE.GetPointY(point))
        Cs_Graph_SPE_3OV.SetPointY(point+32*counter, unsatGraph_SPE.GetPointY(point)/(satGraph_SPE.GetPointY(point)))
        Cs_Graph_SPE_3OV.SetPointError(point+32*counter, satGraph_SPE.GetErrorY(point), unsatGraph_SPE.GetPointY(point)/(satGraph_SPE.GetPointY(point))*np.sqrt((satGraph_SPE.GetErrorY(point)/satGraph_SPE.GetPointY(point))**2+(unsatGraph_SPE.GetErrorY(point)/unsatGraph_SPE.GetPointY(point))**2))
    counter+=1

c = ROOT.TCanvas('c_sodium_calibration','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Na_Graph.SetTitle("Integrated Sodium Charge for Different Digi Settings;Charge at Saturated Setting (pC/keV);Unsaturated/Saturated Integrated Charge")
Na_Graph.SetFillStyle(3001)
Na_Graph.SetMarkerColor(ROOT.kRed)
Na_Graph.SetLineColor(ROOT.kRed)
Na_Graph.GetXaxis().SetLimits(1.0,2.5)
Na_Graph.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph.Draw("ap")
Na_Graph_2OV.SetFillStyle(3001)
Na_Graph_2OV.SetMarkerColor(ROOT.kBlue)
Na_Graph_2OV.SetLineColor(ROOT.kBlue)
Na_Graph_2OV.GetXaxis().SetLimits(1.0,2.5)
Na_Graph_2OV.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph_2OV.Draw("p same")
legend = ROOT.TLegend(0.8, 0.8, 0.9, 0.9)
legend.AddEntry(Na_Graph_2OV, "2 OV")
legend.AddEntry(Na_Graph, "3 OV")
legend.Draw()
c.Print(plotDir+"/"+"sodium_saturation_plot.png")



c = ROOT.TCanvas('c_sodium_calibration_spe','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Na_Graph_SPE.SetTitle("SPE Charge for Different Digi Settings;SPE Charge at Saturated Setting (pC);Unsaturated/Saturated Digi Setting")
Na_Graph_SPE.SetFillStyle(3001)
Na_Graph_SPE.SetMarkerColor(ROOT.kRed)
Na_Graph_SPE.SetLineColor(ROOT.kRed)
Na_Graph_SPE.GetXaxis().SetLimits(2,4)
Na_Graph_SPE.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph_SPE.Draw("ap")
Na_Graph_SPE_2OV.SetFillStyle(3001)
Na_Graph_SPE_2OV.SetMarkerColor(ROOT.kBlue)
Na_Graph_SPE_2OV.SetLineColor(ROOT.kBlue)
Na_Graph_SPE_2OV.GetXaxis().SetLimits(2,4)
Na_Graph_SPE_2OV.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph_SPE_2OV.Draw("p same")
legend = ROOT.TLegend(0.8, 0.8, 0.9, 0.9)
legend.AddEntry(Na_Graph_SPE_2OV, "2 OV")
legend.AddEntry(Na_Graph_SPE, "3 OV")
legend.Draw()
c.Print(plotDir+"/"+"sodium_saturation_spe_plot.png")




c = ROOT.TCanvas('c_sodium_calibration_LO','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Na_Graph_LO.SetTitle("LO for Different Digi Settings;LO Charge at Saturated Setting (pE/MeV);Unsaturated/Saturated Digi Setting")
Na_Graph_LO.SetFillStyle(3001)
Na_Graph_LO.SetMarkerColor(ROOT.kRed)
Na_Graph_LO.SetLineColor(ROOT.kRed)
Na_Graph_LO.GetXaxis().SetLimits(2500,5000)
Na_Graph_LO.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph_LO.Draw("ap")
Na_Graph_LO_2OV.SetFillStyle(3001)
Na_Graph_LO_2OV.SetMarkerColor(ROOT.kBlue)
Na_Graph_LO_2OV.SetLineColor(ROOT.kBlue)
Na_Graph_LO_2OV.GetXaxis().SetLimits(2500, 5000)
Na_Graph_LO_2OV.GetYaxis().SetRangeUser(0.9,1.1)
Na_Graph_LO_2OV.Draw("p same")
legend = ROOT.TLegend(0.8, 0.8, 0.9, 0.9)
legend.AddEntry(Na_Graph_LO_2OV, "2 OV")
legend.AddEntry(Na_Graph_LO, "3 OV")
legend.Draw()
c.Print(plotDir+"/"+"sodium_saturation_LO_plot.png")




linear_fit = ROOT.TF1("Cs_to_Na_SRC_Fit", "[0]*x+[1]")
linear_fit.SetLineColor(ROOT.kBlue)
Cs_Graph_3OV.Fit("Cs_to_Na_SRC_Fit")

c = ROOT.TCanvas('c_sodium_calibration','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Cs_Graph_3OV.SetTitle("Integrated Charge for Different Digi Settings;Cesium Charge at Saturated Setting (pC/keV);Unsaturated Na/Saturated Cs pC/keV Ratio")
Cs_Graph_3OV.SetFillStyle(3001)
Cs_Graph_3OV.SetMarkerColor(ROOT.kRed)
Cs_Graph_3OV.SetLineColor(ROOT.kRed)
Cs_Graph_3OV.GetXaxis().SetLimits(1.0,2.5)
Cs_Graph_3OV.GetYaxis().SetRangeUser(0.9,1.2)
Cs_Graph_3OV.Draw("ap")
#linear_fit.Draw("same")
c.Print(plotDir+"/"+"Cs_to_Na_saturation_plot.png")


c = ROOT.TCanvas('c_sodium_calibration_SPE','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Cs_Graph_SPE_3OV.SetTitle("SPE Charge for Different Digi Settings;SPE Charge at Saturated Setting (pC);SPE Ratio for Different Digi Setting")
Cs_Graph_SPE_3OV.SetFillStyle(3001)
Cs_Graph_SPE_3OV.SetMarkerColor(ROOT.kRed)
Cs_Graph_SPE_3OV.SetLineColor(ROOT.kRed)
Cs_Graph_SPE_3OV.GetXaxis().SetLimits(3.3,4.0)
Cs_Graph_SPE_3OV.GetYaxis().SetRangeUser(0.9,1.2)
Cs_Graph_SPE_3OV.Draw("ap")
c.Print(plotDir+"/"+"Cs_to_Na_saturation_SPE_plot.png")



c = ROOT.TCanvas('c_sodium_calibration_LO','', 1200,1000)
ROOT.gPad.SetGridx()                                                                                                                                                                                                                                                        
ROOT.gPad.SetGridy()
Cs_Graph_LO_3OV.SetTitle("LO for Different Digi Settings;Cs LO at Saturated Setting (pE/MeV);LO Ratio for Different Digi Setting")
Cs_Graph_LO_3OV.SetFillStyle(3001)
Cs_Graph_LO_3OV.SetMarkerColor(ROOT.kRed)
Cs_Graph_LO_3OV.SetLineColor(ROOT.kRed)
Cs_Graph_LO_3OV.GetXaxis().SetLimits(2500,4000)
Cs_Graph_LO_3OV.GetYaxis().SetRangeUser(0.9,1.2)
Cs_Graph_LO_3OV.Draw("ap")
c.Print(plotDir+"/"+"Cs_to_Na_saturation_LO_plot.png")


#code for fitting SRC ratio for Cs -> Na new baseline                                                                                                                                                                                                                           

outfile = ROOT.TFile(plotDir+"/saturation_corrections.root",'RECREATE')
outfile.cd()
linear_fit.Write()
#Na_Graph_2OV.Write()
#Na_Graph.Write()
#Na_Graph_SPE_2OV.Write()
#Na_Graph_SPE.Write()
#Na_Graph_LO_2OV.Write()
#Na_Graph_LO.Write()
#Cs_Graph_3OV.Write()
#Cs_Graph_SPE_3OV.Write()
#Cs_Graph_LO_3OV.Write()
outfile.Close()

