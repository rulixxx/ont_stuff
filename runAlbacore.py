#!/home/sbsuser/.edm/envs/python3/bin/python3

import sys
import os
import re
from random import shuffle
from subprocess import Popen, PIPE
import getopt
import shutil
import h5py
import json
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

vpython  = '/home/sbsuser/.edm/envs/albacore2.3.4/bin/python'
albacore = vpython + ' ~/.edm/envs/albacore2.3.4/bin/'
albacoreVersion = '2.3.4'


def usage() :
   print("""
   runAlbacore.py [-abD] [-j n] [-k kit] [-c conf]  [-f flowcell_type]

      Create chunks of raw fast5 input data to run several instances of 
      the Guppy basecaller in a cluster.

         -j   number of jobs to use (default 40)
         -D   use 1D^2 chemistry
         -b   use the barcode keyword\n\
         -a   get flowcell and kit from fast5 file
         -f   flowcell type (must specify kit too)
         -k   kit used in sequencing (must specify flowcell too)

      Input fast5 must be located in the current directory tree.
      
      Output fastq is located under the stcratch/out directory. A RunInfo json
      is also created with all the experiment's details
""")

try:
   opts, args = getopt.getopt(sys.argv[1:],"k:f:j:abD", ["help"])
except getopt.GetoptError as err:
# print help information and exit
   print(str(err))
   usage()
   exit()

kit = None
flowcellType = None
barcode = False
flowcell = None
start = None
useFast5Info = False
pc = None
runNumber = None
etype = None
minionId = None
minknowV = None
D2 = False
mask = 0o755 #mkdir mask
chunks = 40
pwd = os.getcwd()

for iopt,iarg in opts :
   if iopt in ("-h", "--help"):
      usage()
      sys.exit()
   if iopt == "-j" :
      chunks = int(iarg)
   if iopt == "-b" :
      barcode = True
   if iopt == "-D" :
      D2 = True
   if iopt == "-f" :
      flowcellType = iarg
   if iopt == "-k":
      kit = iarg
   if iopt == "-a" :
      useFast5Info = True

if ( ( kit == None or flowcellType == None) and useFast5Info == False) :
   usage()
   exit()

if ( ( kit != None or flowcellType != None) and useFast5Info == True ) :
   print("cant specify option -a and manually select flowcell and kit")
   usage()
   exit()

if ( ( kit == None and flowcellType == None) and useFast5Info == False ) :
   print("must specify both flowcell and kit")
   usage()
   exit()

#get a list of all the fast5 files 
cmd =  "lfs find %s -type f -name '*.fast5'  2> /dev/null"%pwd
p = Popen(cmd ,shell=True, stdout=PIPE)
(output, error) = p.communicate()
p.wait()
allFast5 = output.decode('utf-8').split('\n')[:-1]

#include any fast5.tmp files in the list
cmd =  "lfs find %s -type f -name '*.fast5.tmp'  2> /dev/null"%pwd
p = Popen(cmd ,shell=True, stdout=PIPE)
(output, error) = p.communicate()
p.wait()
allFast5 += output.decode().split('\n')[:-1]

#read experiment info form a fast5 file
f = h5py.File(allFast5[0],'r')
topGroups = list(f.keys())
if ( 'UniqueGlobalKey' in topGroups ) :
  globalKey = 'UniqueGlobalKey' #single read fast5
else :
  globalKey = list(f.keys())[0] #multi read fast5
flowcellType_fast5 = f['/%s/context_tags/'%globalKey].attrs['flowcell_type'].decode().upper()
kit_fast5 = f['/%s/context_tags/'%globalKey].attrs['sequencing_kit'].decode().upper()
flowcell = f['/%s/tracking_id/'%globalKey].attrs['flow_cell_id'].decode().upper()
start = f['/%s/tracking_id/'%globalKey].attrs['exp_start_time'].decode()
etype = f['/%s/context_tags/'%globalKey].attrs['experiment_type'].decode()
minionId =  f['/%s/tracking_id/'%globalKey].attrs['device_id'].decode()
minknowV =  f['/%s/tracking_id/'%globalKey].attrs['version'].decode()
if (  flowcellType_fast5 == '' or kit_fast5 == '') :
  exit("Couldn't read run information from fast5 file (flowcell type %s, kit %s)!"%(flowcellType_fast5, kit_fast5) )

if ( useFast5Info ) :
   flowcellType= flowcellType_fast5
   kit = kit_fast5
else :
   if ( flowcellType != flowcellType_fast5 or kit != kit_fast5 ) :
      print("Warning specified flowcell %s and kit %s dont match fast5 info ( %s %s"%(flowcellType,kit,flowcellType_fast5,kit_fast5) )

if ( D2 ) :
   chemistry = '1D^2'
else :
   chemistry = '1D'

#store all the run parameters in a json
runInfo = { 'flowcell' : flowcell, 'flowcellType' : flowcellType , 'kit': kit , 'startTime' : start, 'albacoreVersion' : albacoreVersion , 'pc' : pc , 'runNumber' : runNumber, 'barcoded' : barcode, 'experimentType' : etype, 'minionId' : minionId, 'minknowCoreVersion' : minknowV, 'chemistry' : chemistry }
f = open('./RunInfo','w')
json.dump(runInfo,f)
f.close()

#create the scratch directory tree ( if previous run is found remove them)
# ./scratch/out/??/  for outputs
# ./scratch/in/??/   for inputs
for i in range(chunks) :   
   rlabel = "%02d"%i
   rout = "./scratch/out/%s"%(rlabel)
   if ( os.path.exists(rout) ) : shutil.rmtree(rout)
   os.makedirs(rout, mode = mask)
   rin =  "./scratch/in/%s"%(rlabel)
   if ( os.path.exists(rin ) ) : shutil.rmtree(rin)
   os.makedirs(rin, mode = mask)

#create the symlinks
ichunk = 0
cycles = 0
for ifile in allFast5 :
   ifileName = ifile.split('/')[-1]
   ifileName = ifileName.replace('.fast5.tmp','.fast5')
   indir = "./scratch/in/%02d/"%(ichunk)
   src = ifile
   dst = indir + ifileName
   try :
      os.symlink(src, dst)
   except :
      continue
   ichunk += 1   
   if (ichunk == chunks ) :
      cycles += 1
      ichunk = 0

#choose basecaller
if ( D2 == True ) :
   basecaller = 'full_1dsq_basecaller.py'
else :
   basecaller = 'read_fast5_basecaller.py'

#prepare the job array
outfile = open('job.cmd','w')
outfile.write(
"""#!/bin/bash
# @ output = out.%a.log
# @ error = err.%a.log
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 12:59:00
""")
outfile.write("# @ initialdir = .\n")
outfile.write("# @ job_name = %s\n"%(flowcell))
outfile.write("# @ array = 0-%s\n\n"%(chunks-1))
outfile.write("set -e\n")
for ichunk in range(0,chunks)  :
   indir =  "./scratch/in/%02d"%(ichunk)
   outdir = "./scratch/out/%02d"%(ichunk)
   cmd = "if [ $SLURM_ARRAY_TASK_ID -eq %s ]; then\n"%(ichunk )
   cmd +="   export PYTHONPATH='';time %s/%s"%(albacore,basecaller)
   cmd += " -r --flowcell %s --kit %s --output_format fastq --input %s --save_path %s -t 4 --disable_filtering --disable_pings"%(flowcellType,kit,indir,outdir)
   if ( barcode ) : 
      cmd += " --barcoding\n"
   else :
      cmd += "\n"
   cmd += "fi\n"
   outfile.write(cmd)
outfile.close()
launch = '/opt/perf/bin/mnsubmit job.cmd' 
p = Popen(launch, shell=True, stdout=PIPE, stderr=PIPE)
p.wait()
out = p.stdout.read().decode()
err = p.stderr.read().decode()
if ('ERROR' in err) :
   print("!! Error submitting job array to cluster (%s)"%(err))
   raise Exception("Error submitting job array to cluster")
print(out)
