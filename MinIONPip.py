#!/home/sbsuser/.edm/envs/python3/bin/python3

#automatic MinION pipeline for basecalling and demultiplexing
#
#script looks for zip archives in the production directories.
#onece a experiment triggers the script it will:
#
#        * Check experiment with LIMS details
#        * Decompress raw data
#        * Split input reads into batches to allow parallel processing
#        * Submit a job array to do the basecalling with Guppy
#        * Demultiplex fastqs with qcat using a job array
#        * Post stats to the LIMS
#        * Produce QC plots with MinionQC.R
#        * Compress fastq files and remove temporary files

import sys
import os
import re
import time
import random
import json
from subprocess import Popen, PIPE
import requests
import glob
import shutil
import h5py
import pandas as pd
from random import shuffle

#need this otherwise get error when opening fast5 files
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#Divide basecalling amont these many CPUs (about 4hrs to process one flowcell)
chunks = 30

#Guppy basecaller to use
guppyBin= '/home/sbsuser/guppy-cpu-2.2.2/bin/'
guppyVersion = '2.2.2'
#use qcat for demultiplexing
qcat = '/home/sbsuser/.edm/envs/python3/bin/qcat'

#look in these directories for data from the sequencers
minIONdirs = [ '/project/production/shared/minionRAWDATA/dataPC1', '/project/production/shared/minionRAWDATA/dataPC2' ]

#minion run directories
rawData = '/project/production/shared/minionRAWDATA/'

#dir mask
mask = 0o755

#map qcat barcode names to LIMS names
QCAT2LIMS = {}
for i in range(1,25) :
   qbar = "barcode%02i"%(i)
   lbar = "NB%02i"%(i)
   QCAT2LIMS[qbar] = lbar
QCAT2LIMS['none'] = 'unclassified'

#use this pigz
pigz = '/apps/PIGZ/2.3.1/bin/pigz'

#LIMS gmail account
emailUser = 'cnaglims@gmail.com'
emailPwd =  ''
subject = "MinION pipeline error"
warnReceive = [ 'raul.alcantara@cnag.crg.es']
#use this address and headers for interaction with the LIMS
base_url='http://login3'
headers = {'content-type': 'application/json' }
fc_url = "%s%s" % (base_url , "/lims/api/seq/flowcell/" ,  )
lw_url = "%s%s" % (base_url , "/lims/api/seq/loadedwith/" ,  )
fli_url = "%s%s" % (base_url, "/lims/api/seq/flowcell_lane_index" )
subp_url = "%s%s" % (base_url , "/lims/api/seq/subproject/", )

#determine location of the pipeline script so we can import ourselves
whereILive = os.path.dirname(os.path.abspath(sys.argv[0]))

#use this module to send error emails
def send_email( body , subject=subject, recipient=warnReceive ,user=emailUser, pwd=emailPwd ):
   #function that sends email warnings
   import smtplib

   gmail_user = user
   gmail_pwd = pwd
   FROM = user
   TO = recipient if type(recipient) is list else [recipient]
   SUBJECT = subject
   TEXT = body

   # Prepare actual message
   message = """\From: %s\nTo: %s\nSubject: %s\n\n%s
   """ % (FROM, ", ".join(TO), SUBJECT, TEXT)

   try:
      server = smtplib.SMTP("smtp.gmail.com", 587)
      server.ehlo()
      server.starttls()
      server.login(gmail_user, gmail_pwd)
      server.sendmail(FROM, TO, message)
      server.close()
   except:
      raise Exception("failed to send mail!")

#drop state to disk
def updateState(rdir,state) :
   stateFile = rdir + '/state.dat'
   f= open(stateFile,'w')
   json.dump( state, f)
   f.close    
   return


def main() :
   print(minIONdirs)
   for rdir in minIONdirs :
      pc = rdir.split('/')[-1]
      scratchDir = "%s/.%s.scratch"%(rawData,pc)
      stateFile = scratchDir + '/state.dat'

      #look to see if we find any zip archives
      p=Popen("lfs find %s -name '*.zip' | sort -r 2> /dev/null "%(rdir) , shell=True , stdout=PIPE, stderr=PIPE)
      p.wait()
      zipFiles = p.stdout.read().decode().split('\n')[:-1]
      if ( len(zipFiles) == 0 ) : continue
      print("Processing %s"%pc)
      if ( not  os.path.isdir( scratchDir ) ):  
         os.makedirs( scratchDir, mode = mask)
         os.system( "lfs setstripe -c 1 %s"%scratchDir) #change the stripe to 1 for scratch dir
      state = { 'error' : None, 'archives' : {} }
      if ( os.path.isfile(stateFile) ) : 
         f = open(stateFile,'r')
         state = json.load(f)

      #check the size of the zip archives to see if they are finished copying
      #start processing when size of the archive stops changing
      doneCopying = True
      for ifile in zipFiles :
         isize = os.path.getsize( ifile )
         if ( not ifile in state['archives'] ) :
            state['archives'][ifile] = isize
            doneCopying = False
         else :
            if ( state['archives'][ifile] != isize ) :
               state['archives'][ifile] = isize
               doneCopying = False
      if ( not doneCopying ) :
         updateState(scratchDir,state)
         continue

      #get flowcell and kit information one fast5 file
      if (os.path.isfile(rdir+'/RunInfo') ) :
         f= open(rdir+'/RunInfo','r')
         runInfo = json.load(f)
         f.close()
      else :
         p=Popen("unzip -Z1 %s | grep '\.fast5' | head -1 "%(zipFiles[-1]) , shell=True,stdout=PIPE,stderr=PIPE)
         p.wait()
         fast5Path = p.stdout.read().decode().rstrip('\n')
         oneFast5 = "%s/first.fast5"%(scratchDir)
         p=Popen("unzip -p %s %s > %s"%(zipFiles[-1],fast5Path,oneFast5) , shell=True)
         p.wait()
         f = h5py.File(oneFast5,'r')
         topGroups = list(f.keys())
         if ( 'UniqueGlobalKey' in topGroups ) :
            globalKey = 'UniqueGlobalKey' #single read fast5 (old fast5)
         else :
            globalKey = list(f.keys())[0] #multi read fast5 (new fast5)
         flowcellType = f['/%s/context_tags/'%globalKey].attrs['flowcell_type'].decode().upper()
         kit = f['/%s/context_tags/'%globalKey].attrs['sequencing_kit'].decode().upper()
         flowcell = f['/%s/tracking_id/'%globalKey].attrs['flow_cell_id'].decode().upper()
         start = f['/%s/tracking_id/'%globalKey].attrs['exp_start_time'].decode()
         etype = f['/%s/context_tags/'%globalKey].attrs['experiment_type'].decode()
         deviceType =  f['/%s/tracking_id/'%globalKey].attrs['device_type'].decode()
         deviceId =  f['/%s/tracking_id/'%globalKey].attrs['device_id'].decode()
         minknowV =  f['/%s/tracking_id/'%globalKey].attrs['version'].decode()

         runNumber = 1
         chemistry = '1D' #only support 1D chemistry for the time being
         
         os.remove( oneFast5 )
         #dump the run information to a json object
         runInfo = { 'flowcell' : flowcell, 'flowcellType' : flowcellType , 'kit': kit, 'startTime' : start, 'guppyVersion' : guppyVersion , 'pc' : pc , 'runNumber' : runNumber, 'barcoded' : False, 'experimentType' : etype, 'instrumentType' : deviceType, 'instrumentId' : deviceId, 'minknowCoreVersion' : minknowV, 'chemistry' : chemistry }         

         f = open(rdir+'/RunInfo','w')
         json.dump(runInfo,f)
         f.close()

      #raise warning if fast5 doesn't have required information
      if ( runInfo['flowcell'] == '' or runInfo['flowcellType'] == '' or runInfo['kit'] == '') :
         warn01 = "Couldn't read run information from fast5 file at: %s, (flowcell %s, flowcell type %s, kit %s)!"%(runInfo['pc'],runInfo['flowcell'],runInfo['flowcellType'],runInfo['kit'])
         if ( state['error'] != '01' ): 
            send_email(warn01)
            state['error'] = '01'
            updateState(scratchDir,state)
         print(warn01)
         continue

      if ( state['error'] == '01' ) : 
         state['error'] = None
         updateState(scratchDir,state)

      #check cluster queue to see if there are any running jobs for this fc
      mnq = Popen(['/opt/perf/bin/squeue','-o', "%.8i %.40j", "-u", "sbsuser"], stdout=PIPE).communicate()[0].decode()
      if ( "%s_%s_%s"%(pc,runInfo['flowcell'],runInfo['runNumber']) in mnq ) : continue

      #we proceed with decompression once the files have been copyied
      if ( not  os.path.isfile( scratchDir+'/unzip.done') ) :
         print("Unzipping archives from %s"%pc)
         outfile = open(scratchDir + '/unzip.cmd','w')      
         outfile.write(
"""#!/bin/bash
# @ output = unzip.out
# @ error = unzip.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 12:59:00
""")
         outfile.write("# @ initialdir = %s\n"%(scratchDir))
         outfile.write("# @ job_name = %s_%s_%s_unzip\n"%(pc,runInfo['flowcell'],runInfo['runNumber'] ) )
         outfile.write("\n")
         outfile.write("set -e\n")
         for ifile in zipFiles :
            outfile.write("unzip -qq -o %s -d %s\n"%(ifile,scratchDir) )
         outfile.write("touch unzip.done\n")
         outfile.close()
         launch = ['/opt/perf/bin/mnsubmit', scratchDir + '/unzip.cmd' ]
         p = Popen(launch, stdout=PIPE, stderr=PIPE)
         p.wait()
         out = p.stdout.read().decode()
         err = p.stderr.read().decode()
         if ('ERROR' in err) :
            print("!! Error submitting post job to cluster %s (%s)"%(pc,err))
            raise Exception("Error submitting post job to cluster %s"%(pc))
         continue

      #check the error log from the unzip command
      f = open(scratchDir + '/unzip.err','r')
      zipError = f.readlines()
      if ( len(zipError) != 0 ) :
         warn09 = "Zip archive in %s appears corrupted"%(rdir)
         if ( state['error'] != '09' ): 
            state['error'] = '09'
            send_email(warn09)
            updateState(scratchDir,state)
         print(warn09)
         continue
      if ( state['error'] == '09' ) : 
         state['error'] = None
         updateState(scratchDir,state)
    
      #get the LIMS information
      if (os.path.isfile(scratchDir+'/lims.info') ) :
         f= open(scratchDir+'/lims.info','r')
         fliObj = json.load(f)
         f.close()
      else :
         #get the flowcell info from the LIMS
         fc_params = { 'name' :  runInfo['flowcell']  }
         response = requests.get(url=fc_url, headers=headers, params=fc_params)
         if response.status_code != 200 :
            sys.exit("!! Error contacting the LIMS %s"%(response.text))
         fcObj = response.json()['objects']
         if ( len( fcObj ) == 0 ) :
            warn02 = "Flowcell : %s was not found in the LIMS"%(runInfo['flowcell'])
            if ( state['error'] != '02' ) : 
               state['error'] = '02'
               send_email(warn02)
               updateState(scratchDir,state)    
            print(warn02)
            continue
         if ( state['error'] == '02' ) : 
            state['error'] = None
            updateState(scratchDir,state)

         #sort flowcells by run number, assume we are doing the last run
         fcObj.sort(key=lambda x: x['run'], reverse=True)
         runNumber = len( fcObj )
         if (  runInfo['runNumber'] == None ) : 
            runInfo['runNumber'] = runNumber
            f = open(rdir+'/RunInfo','w')
            json.dump(runInfo,f)
            f.close()

         fcStatus = fcObj[0]['flowcell_status']
         #raise an error if the flowcell type doesn't match the LIMS
         if ( fcObj[0]['flowcell_version']['flowcell_version'] != runInfo['flowcellType'] ) :
            warn03 = "Mismatch in flowcell type in dir %s, LIMS: %s , fast5: %s !"%(rdir,fcObj[0]['flowcell_version']['flowcell_version'], runInfo['flowcellType'] )
            if ( state['error'] != '03' ) : 
               state['error'] = '03'
               send_email(warn03)
               updateState(scratchDir,state)
            print(warn03)
            continue
         if ( state['error'] == '03' ) : 
            state['error'] = None
            updateState(scratchDir,state)

         #raise an error if the kit doesn't match the LIMS
         if ( fcObj[0]['kit_version']['kit_version'] != runInfo['kit'] ) :
            warn04 = "Mismatch in kit version in dir %s, LIMS: %s , fast5: %s !"%(rdir,fcObj[0]['kit_version']['kit_version'], runInfo['kit'] )
            if ( state['error'] != '04' ) : 
               state['error'] = '04'
               send_email(warn04)
               updateState(scratchDir,state)  
            print(warn04)
            continue
         if ( state['error'] == '04' ) : 
            state['error'] = None
            updateState(scratchDir,state)

   #      determine if experiment was 1d or 1d2
   #      if ( fcObj[0]['pe_se'] == '/lims/api/seq/pese/6') : d1seq = 1
   #      if ( fcObj[0]['pe_se'] == '/lims/api/seq/pese/7') : d1seq = 2

         #get the FLIs
         fli_params = { 'lane__flowcell__name': runInfo['flowcell'] }
         #fli_parame = { 'lane__flowcell__name': runInfo['flowcell'], 'lane__flowcell__run' : runInfo['runNumber'] } #need to habilitate run for filtering 
         response = requests.get(url=fli_url, headers=headers, params=fli_params)
         if response.status_code != 200 :
            sys.exit("!! Error contacting the LIMS %s"%(response.text))
         fliObj = response.json()['objects']
         nflis = len( fliObj )
        #raise an error if the LIMS doesn't have any FLIs for this experiment
         if ( nflis == 0 ) :
            warn11 = "No FLIs for experiment %s ( %s )"%(runInfo['flowcell'], rdir )
            if ( state['error'] != '11' ) : 
               state['error'] = '11'
               send_email(warn11)
               updateState(scratchDir,state)  
            print(warn11)
            continue
         if ( state['error'] == '11' ) : 
            state['error'] = None
            updateState(scratchDir,state)

         f = open(scratchDir+'/lims.info','w')
         json.dump(fliObj,f)
         f.close()

      #check to see if we have completed all bascalling jobs
      dones = glob.glob("%s/guppy*.done"%scratchDir)
      if ( len(dones) != chunks ) :
         #remove done files from previous execution
         for ifile in dones :
            os.remove( ifile )   
         print("Launching basecalling for %s"%pc)
         #get a list of all the fast5 files 
         cmd =  "lfs find %s -type f -name '*.fast5'  2> /dev/null "%scratchDir
         p = Popen(cmd ,shell=True, stdout=PIPE)
         (output, error) = p.communicate()
         p.wait()
         allFast5 = output.decode().split('\n')[:-1]

         #include any fast5.tmp files in the list
         cmd =  "lfs find %s -type f -name '*.fast5.tmp'  2> /dev/null "%scratchDir
         p = Popen(cmd ,shell=True, stdout=PIPE)
         (output, error) = p.communicate()
         p.wait()
         allFast5 += output.decode().split('\n')[:-1]

         shuffle( allFast5 ) #to evenout execution times of all the jobs


         #create the directory tree, if previous run is found remove it
         for i in range(chunks) :   
            rlabel = "%02d"%i
            rout = scratchDir+ "/scratch/out/%s"%(rlabel)
            if ( os.path.exists(rout) ) : shutil.rmtree(rout)
            os.makedirs(rout, mode = mask)
            rin = scratchDir+ "/scratch/in/%s"%(rlabel)
            if ( os.path.exists(rin ) ) : shutil.rmtree(rin)
            os.makedirs(rin, mode = mask)

         #create the symlinks
         ichunk = 0
         cycles = 0
         for ifile in allFast5 :
            ifileName = ifile.split('/')[-1]
            ifileName = ifileName.replace('.fast5.tmp','.fast5')
            indir = "%s/scratch/in/%02d/"%(scratchDir,ichunk)
            src = ifile
            dst = indir + ifileName
            os.symlink(src, dst)
            ichunk += 1   
            if (ichunk == chunks ) :
               cycles += 1
               ichunk = 0
        
         if ( len( fliObj ) > 1 ) : 
            runInfo['barcoded'] = True
            f = open(rdir+'/RunInfo','w')
            json.dump(runInfo,f)
            f.close()
         
         #prepare the job array
         outfile = open(scratchDir + '/job.cmd','w')
         outfile.write(
"""#!/bin/bash
# @ output = guppy.%a.out
# @ error = guppy.%a.err
# @ total_tasks = 1
# @ cpus_per_task = 8
# @ wall_clock_limit = 24:59:00
""")
         outfile.write("# @ initialdir = %s\n"%(scratchDir))
         outfile.write("# @ job_name = %s_%s_%s_base\n"%(pc,runInfo['flowcell'],runInfo['runNumber']))
         outfile.write("# @ array = 0-%s\n\n"%(chunks-1))
         #load the required libraries
         outfile.write("module purge\n")
         outfile.write("module load gcc/4.9.3-gold\n")
         outfile.write("module add hdf5/1.8.12\n")
         outfile.write("export LD_LIBRARY_PATH=/home/sbsuser/guppy-cpu-2.2.2/lib/:$LD_LIBRARY_PATH\n")
         outfile.write("set -e\n")
         for ichunk in range(0,chunks)  :
            indir =  "./scratch/in/%02d"%(ichunk)
            outdir = "./scratch/out/%02d"%(ichunk)
            cmd = "if [ $SLURM_ARRAY_TASK_ID -eq %s ]; then\n"%(ichunk )
            cmd +="   time %s/guppy_basecaller"%(guppyBin)
            cmd +="  --flowcell %s --kit %s"%(runInfo['flowcellType'],runInfo['kit'])
            cmd += " --input %s --save_path %s --runners 8 -t 1 --disable_pings -r\n"%(indir,outdir)
            cmd += "  touch guppy.%d.done\n"%(ichunk)
            cmd += "fi\n\n"
            outfile.write(cmd)
         outfile.close()
         launch = ['/opt/perf/bin/mnsubmit', scratchDir + '/job.cmd' ]
         p = Popen(launch, stdout=PIPE, stderr=PIPE)
         p.wait()
         out = p.stdout.read().decode()
         err = p.stderr.read().decode()
         if ('ERROR' in err) :
            print("!! Error submitting job array to cluster %s (%s)"%(pc,err))
            raise Exception("Error submitting job array to cluster %s"%(pc))
            return
         continue  

      #proceed to demultiplex if experiment is barcoded
      fastqDir = scratchDir +"/scratch/out"
      if ( runInfo['barcoded'] ) :
         dones = glob.glob("%s/demux*.done"%scratchDir)
         if ( len(dones) != chunks ) :
            print("Launching demux for %s"%pc)
            cmd =  "lfs find %s -type f -name '*.fastq'  2> /dev/null"%(fastqDir)
            p = Popen(cmd ,shell=True, stdout=PIPE)
            (output, error) = p.communicate()
            p.wait()
            allFastq = output.decode().split('\n')[:-1]
            if ( len( allFastq ) == 0 ) : sys.exit("Didn't find any fastq files in %s"%fastqDir)

            #try to even out execution times
            shuffle( allFastq )
            #create the demux scratch directory tree 
            # ./scratch/demux/in/??  for input fastqs
            # ./scratch/demux/out/?? for summary outputs
            for i in range(chunks) :   
               rlabel = "%02d"%i
               rout = "%s/scratch/demux/out/%s"%(scratchDir,rlabel)
               if ( os.path.exists(rout) ) : shutil.rmtree(rout)
               os.makedirs(rout, mode = mask)
               rin =  "%s/scratch/demux/in/%s"%(scratchDir, rlabel)
               if ( os.path.exists(rin ) ) : shutil.rmtree(rin)
               os.makedirs(rin, mode = mask)

            #create symlinks to the input fastqs
            ichunk = 0
            cycles = 0
            for ifile in allFastq :
               ifileName = ifile.split('/')[-1]
               indir = "%s/scratch/demux/in/%02d/"%(scratchDir,ichunk)
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

            #set the demultiplexing qit to use with qcat
            if ( len( fliObj ) == 12 ) :
               demuxKit = 'NBD103/NBD104' #12 barcodes
            else :
               demuxKit = 'NBD104/NBD114' #24 barcodes

            #prepare the job array for the cluster
            outfile = open(scratchDir + '/demux.cmd','w')
            outfile.write(
"""#!/bin/bash
# @ output = demux.%a.out
# @ error = demux.%a.err
# @ total_tasks = 1
# @ cpus_per_task = 1
# @ wall_clock_limit = 0:30:00
""")
            outfile.write("# @ initialdir = %s\n"%(scratchDir))
            outfile.write("# @ job_name = %s_%s_%s_demux\n"%(pc,runInfo['flowcell'],runInfo['runNumber']))
            outfile.write("# @ array = 0-%s\n\n"%(chunks-1))
            #load the required libraries
            outfile.write("export PYTHONPATH=''\n")
            outfile.write("set -e\n")
            for ichunk in range(0,chunks)  :
               indir =  "./scratch/demux/in/%02d"%(ichunk)
               outdir = "./scratch/demux/out/%02d"%(ichunk)
               cmd = "if [ $SLURM_ARRAY_TASK_ID -eq %s ]; then\n"%(ichunk )
               cmd += "   cat %s/*.fastq | %s -k %s -b %s\n"% (indir,qcat,demuxKit,outdir)
               cmd += "   touch demux.%d.done\n"%(ichunk)
               cmd += "fi\n\n"
               outfile.write(cmd)
            outfile.close()
            launch = ['/opt/perf/bin/mnsubmit', scratchDir + '/demux.cmd' ] 
            p = Popen(launch, stdout=PIPE, stderr=PIPE)
            p.wait()
            out = p.stdout.read().decode()
            err = p.stderr.read().decode()
            if ('ERROR' in err) :
               print("!! Error submitting job array to cluster %s (%s)"%(pc,err))
               raise Exception("Error submitting job array to cluster %s"%(pc))
               return
            continue

      #post processing. accumulate fastq and summary files, split summary files and run minionqc script 
      if ( not  os.path.isfile( scratchDir+'/post.done') ) :

         print("Launching postprocessing for %s"%pc)
         #check the barcodes
         if ( runInfo['barcoded'] ) :
            cmd =  "lfs find %s/scratch/demux -name 'barcode*.fastq'  2> /dev/null "%(scratchDir)
            p = Popen(cmd ,shell=True, stdout=PIPE)
            (output, error) = p.communicate()
            p.wait()
            fout = output.decode().split('\n')[:-1]         
            barcodes = set()
            for iout in fout :
               ibar = iout.split('/')[-1].split('.')[0]
               barcodes.add( ibar )
         #setup script 
         outfile = open(scratchDir + '/post.cmd','w')
         outfile.write(
"""#!/bin/bash
# @ output = post.out
# @ error = post.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 03:59:00
""")
         outfile.write("# @ initialdir = %s\n"%(scratchDir))
         outfile.write("# @ job_name = %s_%s_%s_post\n"%(pc,runInfo['flowcell'],runInfo['runNumber'] ) )
         outfile.write("\n")
         outfile.write("module purge\n")
         outfile.write("module load gcc/6.3.0\n")
         outfile.write("module load parallel/20150222\n")
         outfile.write("module load R/3.5.0\n")
         outfile.write("export R_LIBS_USER=\"/home/sbsuser/MinION/R\"\n\n")
         outfile.write("set -e\n\n")

         #merge all the fastq files
         if (  runInfo['barcoded'] ) :
            for ibarcode in ( list(barcodes) + ['none'] )  :
               fliName = "%s_%s_%s"%(runInfo['flowcell'],runInfo['runNumber'], QCAT2LIMS[ibarcode])
               fastqFile = "%s/%s_%s_%s.fastq"%(rdir,runInfo['flowcell'],runInfo['runNumber'], QCAT2LIMS[ibarcode] )            
               outfile.write("lfs find %s -name '*.fastq'  | grep %s | xargs -i cat {} > %s\n"%(scratchDir + "/scratch/demux/out",ibarcode,fastqFile) )
         else :
            fliName = "%s_%s_%s"%(runInfo['flowcell'],runInfo['runNumber'], '0')
            fastqFile = "%s/%s_%s_0.fastq"%(rdir,runInfo['flowcell'],runInfo['runNumber'])
            outfile.write("lfs find %s -name '*.fastq' | xargs -i cat {} > %s\n"%(scratchDir + "/scratch/out",fastqFile) )

         #merge all the sequencing summary files
         outfile.write("lfs find %s -name 'sequencing_summary.txt' | head -1 | xargs head -1 > %s/sequencing_summary.0.txt\n"%(scratchDir +'/scratch/out',rdir))
         outfile.write("lfs find %s -name 'sequencing_summary.txt' | xargs -i tail -n+2 {} >> %s/sequencing_summary.0.txt\n"%(scratchDir +'/scratch/out',rdir))
         outfile.write("cd %s\n"%rdir)
         if (  runInfo['barcoded'] ) :
            #split summary_files by barcode
            #a bit messy but this manages to run parallel instances of the filter function
            outfile.write("ls *.fastq | parallel -j4 \" /home/sbsuser/.edm/envs/python3/bin/python3 -c \\\"import sys;sys.path.append('%s');import MinIONPip;MinIONPip.filterSS('{}')\\\" \"\n"%(whereILive))
            outfile.write("rm sequencing_summary.0.txt\n")
            for ibarcode in ( list(barcodes) )  :
               #generate MinionQC plots
               fliName = "%s_%s_%s"%(runInfo['flowcell'],runInfo['runNumber'], QCAT2LIMS[ibarcode])
            outfile.write("ls *.fastq | cut -d'_' -f 3 | cut -d'.' -f 1 | parallel -j 4 \"Rscript /home/sbsuser/MinION/minion_qc-1.3.0/MinionQC.R -f %s_{} -i %s -o %s/qc_plots/{}\"\n"%(runInfo['flowcell'],'./sequencing_summary.{}.txt','.'))
         else :
            outfile.write("Rscript /home/sbsuser/MinION/minion_qc-1.3.0/MinionQC.R -f %s -i %s -o %s/qc_plots/0\n"%(fliName,'./sequencing_summary.0.txt',scratchDir))

         #touch post done to show we completed successfully 
         outfile.write("touch %s/post.done\n\n"%scratchDir)
         outfile.close()
         launch = ['/opt/perf/bin/mnsubmit', scratchDir + '/post.cmd' ] 
         p = Popen(launch, stdout=PIPE, stderr=PIPE)
         p.wait()
         out = p.stdout.read().decode()
         err = p.stderr.read().decode()
         if ('ERROR' in err) :
            print("!! Error submitting post job to cluster %s"%(pc))
            raise Exception("Error submitting post job to cluster %s"%(pc))
         continue

      #consistency check 1
      summaryFiles = glob.glob("%s/sequencing_summary.*.txt"%rdir)
      fastqFiles = glob.glob("%s/*.fastq"%rdir)
      #raise warning if number of fastq and summary files dont match
      if ( len(summaryFiles) != len(fastqFiles) ) :
         warn06 = "Disagreeing number of summary and fastq files in %s ( %s )!"%(runInfo['flowcell'],runInfo['pc'])
         if ( state['error'] != '06' ): 
            send_email(warn06)
            state['error'] = '06'
            updateState(scratchDir,state)
         print(warn06)
         continue
      if ( state['error'] == '06' ) : 
         state['error'] = None
         updateState(scratchDir,state)

      #get the number of reads and yield from the fastq (might differ from the numbers obtained from sequencing_summary)
      # for consistency check 2
      tot_yyield = 0
      tot_nreads = 0
      if ( not runInfo['barcoded'] ) :
         QCAT2LIMS['0'] = '0'
         barcodes = {'0'}
      cmd = "cat %s/*.fastq | wc -c -l"%rdir
      p=Popen( cmd, shell=True , stdout=PIPE, stderr=PIPE)
      p.wait()
      out = p.stdout.read().decode()
      err = p.stderr.read().decode()
      out = out.rstrip('\n')
      tot_nreads = int(out.split(' ')[-2])//4
      tot_yyield = int(out.split(' ')[-1])//4 - tot_nreads
      
      print("Posting data to the LIMS")
      tot_nreads2 = 0 #compute number again from the summary files
      stats = {} #store yield and nreads and other stats
      for isummaryFile in summaryFiles  :
         ibarcode = isummaryFile.split('.')[1]
         oo = {'mean_qscore_template' : [] , 'sequence_length_template' : [] }
         f = open(isummaryFile,'r')
         head = f.readline()
         imq = head.split().index( 'mean_qscore_template' )
         islt = head.split().index( 'sequence_length_template' )
         for iline in f :
            oo['mean_qscore_template'].append( float(iline.split('\t')[imq] ))
            oo['sequence_length_template'].append( int(iline.split('\t')[islt]) )
         f.close()
         #compute n. reads; yield; N50; max median mean fragment size
         df = pd.DataFrame(oo)
         df = df.sort_values(by='sequence_length_template')
         yyield = df['sequence_length_template'].sum()
         nreads = len( df )
         tot_nreads2 += nreads
         cumsum = 0
         for i in df['sequence_length_template'] :
             cumsum += i
             if ( cumsum > yyield/2.0) : break
         N50 = i
         medians = df.median()
         means = df.mean()
         maxs = df.max()

         #post yield and number of clusters to LIMS
         ncounts = str( int(round( nreads /1000.0) ) )
         pct_total_reads = nreads * 100.0/tot_nreads
         yieldpfStr = str( int(round(yyield/1.0e6) ))
         lw_payload = { 'clusters_pass_filt' : ncounts, "percent_total_reads": pct_total_reads, "yield_pass_filt": yieldpfStr }
         #find the object in fliObj that has our index
         if ( ibarcode == '0' ) :
            ifliObj =  fliObj[0] 
         else :
            for ifliObj in fliObj :
               if ( ifliObj['index_name'] == ibarcode ) : break 
            else :
               continue

         if ( ibarcode != 'unclassified' ) :
            print("%s %s"%(isummaryFile, lw_payload))
            url = "%s/%s" % (lw_url , ifliObj["id"],  )
            response = requests.patch(url=url,  headers=headers, data=json.dumps(lw_payload))
            if response.status_code != 202 :
               print("!! Error in patching the LIMS %s"%(response.text))
               raise Exception(" !! Error in patching the LIMS %s"%(response.text) )

            #save stats by lwid to produce a report
            stats[ ifliObj["id"] ] = { 'yield' : yyield, 'reads' : nreads , 'N50' : N50, 'max' : maxs, 'medians' : medians, 'means' : means }
         else :
            stats[ 'unclassified' ] = { 'yield' : yyield, 'reads' : nreads , 'N50' : N50, 'max' : maxs, 'medians' : medians, 'means' : means }
         
         if ( ibarcode == '0' ) : break #if experiment is not barcoded just loop once

      #check fastq files against summary files consistency check 2
      if ( tot_nreads != tot_nreads2 ) : 
         warn10 = "Mismatch between reads in fastq files (%s) and summary files (%s) at  %s "%(tot_nreads,tot_nreads2,rdir)
         if ( state['error'] != '10' ): 
            state['error'] = '10'
            send_email(warn10)
            updateState(scratchDir,state)
         print(warn10)
         continue
      if ( state['error'] == '10' ) : 
         state['error'] = None
         updateState(scratchDir,state)

      #copy qc plots to LIMS accesible directory
      qcDirTarget = '/project/lims/minion_qc/%s/%s/'%(runInfo['flowcell'],runInfo['runNumber'])
      qcDirLink = "http://lims.cnag.cat/public/minion_qc/%s/%s/"%(runInfo['flowcell'],runInfo['runNumber'])
      if ( not os.path.isdir( qcDirTarget ) ): os.makedirs( qcDirTarget, mode = mask )
      p=Popen("cp -r %s %s" % ( rdir+'/qc_plots/*', qcDirTarget ), shell=True)
      p.wait()

      #create final report and final directory name
      f = open(rdir+'/%s_%s_summary.csv'%(runInfo['flowcell'],runInfo['runNumber']),'w')
      f.write('\n')
      f.write('Sample details\n')
      f.write('subproject,sample_name,sample_barcode,aliquot,library,index,reads,yield_[bp],N50,max_length,mean_length,median_length,mean_q,median_q\n')
      try :
         expDate = runInfo['startTime'].split('T')[0].replace('-','')
      except :
         expDate = fcObj[0]['date'].replace('-','')

      for ifliObj in fliObj :
         subp = ifliObj['subproject_name']
         libn = ifliObj['library_name']
         libraryURL = "%s%s" % (base_url , ifliObj['library'] ,  )
         response = requests.get(url=libraryURL , headers=headers,)
         if response.status_code != 200 :
            sys.exit("Problem contacting LIMS server" , libraryURL)
         libObj = response.json()
         aliquotURL = "%s%s" % (base_url , libObj['aliquot'] ,  )
         response = requests.get(url=aliquotURL , headers=headers,)
         if response.status_code != 200 :
            sys.exit("Problem contacting LIMS server" , sampleURL)
         aliquotObj = response.json()
         sampleURL = "%s%s" % (base_url , libObj['sample'] ,  )
         response = requests.get(url=sampleURL , headers=headers,)
         if response.status_code != 200 :
            sys.exit("Problem contacting LIMS server" , sampleURL)
         sampleObj = response.json()
         aliquot = aliquotObj['aliquot_barcode']
         sampb = sampleObj['barcode']
         sampn = sampleObj['name']
         istats = stats[ ifliObj["id"] ]
         f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(subp,sampn,sampb,aliquot,libn,ifliObj['index_name'],istats['reads'], istats['yield'], istats['N50'],istats['max']['sequence_length_template'],istats['means']['sequence_length_template'],istats['medians']['sequence_length_template'],istats['means']['mean_qscore_template'],istats['medians']['mean_qscore_template'] )) 
      if ( 'unclassified' in stats ) :
          f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%('','','','','','unclassified',istats['reads'], istats['yield'], istats['N50'],istats['max']['sequence_length_template'],istats['means']['sequence_length_template'],istats['medians']['sequence_length_template'],istats['means']['mean_qscore_template'],istats['medians']['mean_qscore_template'] ))
      f.write('\n')
      f.write('\n')
      f.write('Experiment details\n')
      f.write('flowcell,%s\n'%runInfo['flowcell'])
      f.write('flowcell_type,%s\n'%runInfo['flowcellType'])
      f.write('kit,%s\n'%runInfo['kit'])
      f.write('experiment_type,%s\n'%runInfo['experimentType'])
      f.write('MinKnow_version,%s\n'%runInfo['minknowCoreVersion'])
      f.write('Device,%s\n'%runInfo['instrumentType'])
      f.write('Guppy_version,%s\n'%runInfo['guppyVersion'])
      f.close()

      if ( runInfo['barcoded'] ) :
         dirName = "%s_%s_barcoded_%s_%s_%s"%(expDate,subp,pc.replace('data',''),runInfo['runNumber'],runInfo['flowcell'])
      else :
         #experiment name, subproject, library name, sample barcode, pc id, run Numeber and flowcell
         dirName = "%s_%s_%s_%s_%s_%s_%s"%(expDate,subp,libn,sampb,pc.replace('data',''),runInfo['runNumber'],runInfo['flowcell'])

      if ( dirName == '' ) :
         warn08 = "Bad directory name %s"%(rdir)
         if ( state['error'] != '08' ): 
            state['error'] = '08'
            send_email(warn08)
            updateState(scratchDir,state) 
         print(warn08)
         continue 
      if ( state['error'] == '08' ) : 
         state['error'] = None
         updateState(scratchDir,state)

      #try to move files to final directory
      targetDir = rawData + dirName
      if ( os.path.isdir( targetDir ) ):
         warn07 = "Directory %s already exists! Can't move run from %s"%(targetDir, pc)
         if ( state['error'] != '07' ): 
            state['error'] = '07'
            send_email(warn07)
            updateState(scratchDir,state) 
         print(warn07)
         continue 
      else :
          os.makedirs( targetDir, mode = mask )
      if ( state['error'] == '07' ) : 
         state['error'] = None
         updateState(scratchDir,state)
      p=Popen("mv %s %s" % ( rdir+'/*', targetDir), shell=True)
      p.wait()
      #create the clean directory and move all temp files there
      #this way we can start processing a new run while still doing file cleanup
      cleanDir = scratchDir + '.clean'
      if ( os.path.isdir( cleanDir ) ): shutil.rmtree( cleanDir )
      os.makedirs( cleanDir, mode = mask )
      p=Popen("mv %s %s" % ( scratchDir+'/*', cleanDir), shell=True)
      p.wait()
      print("Launching cleaning script for %s"%pc)
      outfile = open(cleanDir + '/clean.cmd','w')
      outfile.write(
"""#!/bin/bash
# @ output = clean.out
# @ error = clean.err
# @ total_tasks = 1
# @ cpus_per_task = 8
# @ wall_clock_limit = 5:59:00
""")
      outfile.write("# @ initialdir = %s\n"%(cleanDir))
      outfile.write("# @ job_name = %s_%s_%s_clean\n"%(pc,runInfo['flowcell'],runInfo['runNumber']))
      outfile.write("\n")
      outfile.write("module load pigz\n")
      outfile.write("module load parallel\n")
      outfile.write("set -e\n\n")
      outfile.write("pigz -p 8 %s/*.fastq\n"%targetDir)
      outfile.write("lfs find %s -type f | grep -v 'clean\.' | parallel -j3 rm\n"%(cleanDir ))
      outfile.write("lfs find %s -type l |  parallel -j3 rm\n"%(cleanDir ))
      outfile.write("lfs find %s | grep -v 'clean\.' | sort -r -u | xargs rm -rf\n"%(cleanDir ))
      outfile.write("touch clean.done\n")
      outfile.close()
      launch = ['/opt/perf/bin/mnsubmit', cleanDir + '/clean.cmd' ]
      p = Popen(launch, stdout=PIPE, stderr=PIPE)
      p.wait()
      out = p.stdout.read().decode()
      err = p.stderr.read().decode()
      if ('ERROR' in err) :
         print("!! Error submitting clean job to cluster %s"%(pc))
         raise Exception("Error submitting clean job to cluster %s"%(pc))
      continue

#call this on cluster to split the sequencing summary file multiplexed runs
# being that qcat doenst do this splitting on its own
def filterSS( fastq ) :
   try :
      ff = open( fastq, 'r')
   except FileNotFoundError:
      raise Exception( ' File %s not found !'%fastq )

   try :
      sf = open( 'sequencing_summary.0.txt' , 'r')
   except FileNotFoundError:
      raise Exception( ' sequencing_summary.0.txt not found !' )

   matchOb = re.search('_([^_]+)\.fastq',fastq )
   indexName = matchOb.group(1)
   outFile = 'sequencing_summary.%s.txt'%indexName
   nline = 0
   ids = set()
   for iline in ff :
      nline += 1
      if ( nline % 4 != 1 ) : continue
      ids.add( iline.split(' ')[0][1:] )
   ff.close()

   of = open(outFile,'w')
   head = sf.readline()
   of.write(head)
   for iline in sf :
      iid = iline.split('\t')[1]
      if ( iid in ids ) : of.write(iline)
   of.close()


if __name__ == "__main__":
   main()
      
