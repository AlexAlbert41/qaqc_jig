#! /usr/bin/python3

import os
import sys
import time
import numpy as np
from optparse import OptionParser
from datetime import datetime
from simple_pid import PID
from subprocess import Popen, PIPE
import subprocess
import logging
from Chillers import *
from Temperatures import *


parser = OptionParser()

parser.add_option("--target", type=float, dest="target", default=23.0)
parser.add_option("--initialDelay", type=float, dest="initialDelay", default=30)
parser.add_option("--initialTemp", type=float, dest="initialTemp", default=21.0)
parser.add_option("--debug", action="store_true")
(options, args) = parser.parse_args()

debug = False
if options.debug:
    debug = True


now = datetime.now()
this_time = now.strftime('%Y-%m-%d_%H:%M:%S')
logfile = os.path.expanduser('~/.setBoxTemp_PID_%s.log'%this_time)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S',filename=logfile,level=logging.INFO)

min_temp = options.target - 6.
max_temp = options.target + 6. #interval needs to be symmetrical not to diverge when a temp spike occurs

min_temp_safe = 10.
max_temp_safe = 40.

if options.target < min_temp_safe or options.target > max_temp_safe:
    print("### ERROR: set temperature outside allowed range ["+str(min_temp_safe)+"-"+str(max_temp_safe)+"]. Exiting...")
    sys.exit(-1)


mychiller = Chiller()
state = mychiller.getState()

logging.info(">>> Chiller::state: "+str(state))

if int(state) == 0:
    logging.info("--- powering on the chiller")
    mychiller.setState('1')
    time.sleep(5)
    state = mychiller.getState()
    logging.info(">>> Chiller::state: "+str(state))
    if state == 0:
        logging.info("### ERROR: chiller did not power on. Exiting...")
        sys.exit(-2)


box_temp = 23.
while True:
    try:
        box_temp = np.mean(read_box_temp())
        break
    except ValueError as e:
        logging.info(e)
        time.sleep(5)
        continue


water_temp = mychiller.getWaterTemp()
new_temp = options.initialTemp
logging.info("--- setting chiller water temperature at "+str(round(float(new_temp),2))+" C   [box temperature: "+str(round(float(box_temp),2))+" C   water temperature: "+str(round(float(water_temp),2))+" C]")

mychiller.setTemp(new_temp)
sleep_time = options.initialDelay
logging.info("--- sleeping for "+str(sleep_time)+" s\n")
sys.stdout.flush()
time.sleep(sleep_time)


pid = PID(0.2, 0., 50., setpoint=options.target)
pid.output_limits = (min_temp-options.target, max_temp-options.target)



while True:
    try:
        try:
            box_temp = np.mean(read_box_temp())
        except ValueError as e:
            logging.info(e)
            time.sleep(5)
            continue
        
        output = pid(box_temp)
        new_temp += output
        
        #safety check
        new_temp = min([max([new_temp,min_temp]),max_temp])
        
        if debug:
            p, i, d = pid.components
            logging.info("== DEBUG == P=", p, "I=", i, "D=", d)
        
        water_temp = mychiller.getWaterTemp()
        logging.info("--- setting chiller water temperature at "+str(round(float(new_temp), 2))+" C   [box temperature: "+str(round(float(box_temp),2))+" C   water temperature: "+str(round(float(water_temp),2))+" C]")

        mychiller.setTemp(new_temp)
        sleep_time = 31
        logging.info("--- sleeping for "+str(sleep_time)+" s   [kill at any time with ctrl-C]\n")
        sys.stdout.flush()
        time.sleep(sleep_time)
        
    except KeyboardInterrupt:
        break

logging.info("--- powering off the chiller")
mychiller.setState('0')
time.sleep(5)
state = mychiller.getState()
logging.info(">>> Chiller::state: "+str(state))
if state == 0:
    logging.info("### ERROR: chiller did not power off. Exiting...")
    sys.exit(-3)
logging.info("bye")
sys.exit(0)
