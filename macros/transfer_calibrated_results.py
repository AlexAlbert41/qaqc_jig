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
 
 
input_path = '/data/QAQC_SM/qaqc-gui_output/Sodium_Production_ReCalibration/'
output_path = '/run/media/cptlab/qaqc_data/results/QAQC_SM/qaqc-gui_output/SM_QAQC_Production'
#plotDir = '/data/QAQC_SM/qaqc-gui_output/SummaryPlots_Calibrations_first165/'
gui_settings_path = '/run/media/cptlab/qaqc_data/results/QAQC_SM/qaqc-gui_output/SM_QAQC_Production'

files_to_transfer = glob.glob(input_path+'/*_analysis_new_sodium_calib.root')
uncalib_files_paths = glob.glob(output_path+'/run*/*_analysis.root')
print(files_to_transfer)
#print(uncalib_files_paths)

for input_file in files_to_transfer:
    tokens = input_file.split('/')
    module = ''
    for token in tokens:
        if 'module' in token:
            module = token[6:20]
    #print(module)
    for uncalib_root in reversed(uncalib_files_paths): #reversed because the generated calibrated root file was the from the most recent data taking run for that module
        if module not in uncalib_root:
            continue
        copy_dir = uncalib_root[:uncalib_root.rfind('/')+1]
        #print(f"cp {input_file} {copy_dir}/module_{module}_analysis_both_calibs_satCorrection.root")
        os.system(f"cp {input_file} {copy_dir}/module_{module}_analysis_recalibration.root")
        integrals_dir = input_file.replace(".root","")
        integrals_dir = integrals_dir.replace("analysis", "integrals")
        integrals_dir_end  = integrals_dir[integrals_dir.rfind('/')+1:]
        integrals_dir_end = integrals_dir_end.replace("module", "module_")
        #print(f"cp -r {integrals_dir} {copy_dir}/{integrals_dir_end}")
        os.system(f"cp -r {integrals_dir} {copy_dir}/{integrals_dir_end}")
