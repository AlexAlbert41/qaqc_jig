

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

#data_path = '/home/cptlab/mnt/btl-upload/upload/QAQC_SM/qaqc-gui_output/SM_QAQC_Production/'
data_path = '/home/cptlab/mnt/btl-upload/upload/QAQC_SM/qaqc-gui_output/SM_QAQC_Production'
data_path_2 = '/home/cptlab/mnt/btl-upload/upload/from_cptlab_100724/qaqc-gui_output/SM_QAQC_Production/'
selections = []
plotDir = '/data/QAQC_SM/qaqc-gui_output/Sodium_Production_ReCalibration/'


#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.25,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)



modules = []
modules_int = []
params = {}
# retrieving root files 
inputFiles1 = glob.glob(data_path+'/run*/*_integrals.hdf5')
gui_settings_files1 = glob.glob(data_path+'/run*/qaqc_gui.settings')
#inputFiles2 = glob.glob(data_path_2+'/run*/*_integrals.hdf5')
gui_settings_files2 = glob.glob(data_path+'/run*/qaqc_gui.settings')
#inputFiles = inputFiles1+inputFiles2
inputFiles = inputFiles1
gui_settings_files = gui_settings_files1+gui_settings_files2
print(inputFiles)
print(gui_settings_files)
for inputFile in inputFiles:
    tokens = inputFile.split('/')
    run = ''
    for token in tokens:
        if 'module' in token:
            module = token[7:21] # SM ID
        if 'run' in token:
            run = int(token[3:]) # run number
    modules.append(module)
    modules_int.append(int(module))
    params[module] = [inputFile,run,'GOOD']
modules_int.sort()
print(modules_int)

if not os.path.isdir(plotDir):
    os.mkdir(plotDir)

print(len(modules))

modules_tested = ["32110020000016"]
counter=0
for num, module in enumerate(modules):
    #print("Module: ", module)
    #print(counter)
    counter+=1
    #if counter>1:
    #    break
    if module in modules_tested:
        continue
    modules_tested.append(module)
    param = params[module]
    accept = 1
    #if "32110020008617" not in module:continue
    #print(int(module))
    if int(module)<32110020008694:continue #this is when we started using sodium 
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
    #{plotDir}    accept *= tempAccept
    if accept == 0:
        continue
    out_file = f"{plotDir}/module{module}_analysis_new_sodium_calib.root"
    #print(f"~/AlexAlbert41/qaqc_jig/python/analyze-waveforms {params[module][0]} -o {out_file} --sourceType cesium --print-pdfs {plotDir}")
    os.system(f"analyze-waveforms {params[module][0]} -o {out_file} --slot {slot} --calibration 'both' --sourceType sodium --print-pdfs {plotDir}")
