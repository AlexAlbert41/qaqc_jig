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
 
 
input_path = '/data/QAQC_SM/qaqc-gui_output/SM_results_after_calibrations_scale_and_channel_copy/'
output_path = '/data/QAQC_SM/qaqc-gui_output/SM_QAQC_Production'
plotDir = '/data/QAQC_SM/qaqc-gui_output/SummaryPlots_Calibrations_first165/'
gui_settings_path = '/data/QAQC_SM/qaqc-gui_output/SM_QAQC_Production/'

files_to_transfer = glob.glob(input_path+'/*_analysis__both_calibs.root')
uncalib_files_paths = glob.glob(output_path+'/run*/*_analysis.root')

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
        #print(f"cp {input_file} {copy_dir}/module_{module}_analysis_both_calibs.root")
        os.system(f"cp {input_file} {copy_dir}/module_{module}_analysis_both_calibs.root")
        integrals_dir = input_file.replace("module", "module_")
        integrals_dir = integrals_dir.replace("analysis__both_calibs.root", "integrals")
        integrals_dir_end  = integrals_dir[integrals_dir.rfind('/')+1:]
        #print(f"cp -r {integrals_dir} {copy_dir}/{integrals_dir_end}_both_calibs")
        os.system(f"cp -r {integrals_dir} {copy_dir}/{integrals_dir_end}_both_calibs")
