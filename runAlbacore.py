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
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE" #caused error w/o this

vpython  = '/home/sbsuser/.edm/envs/albacore2.3.4/bin/python'
albacore = vpython + ' /home/sbsuser/.edm/envs/albacore2.3.4/bin/'
albacoreVersion = '2.3.4'
pwd = os.getcwd()
mask = 0o755 #mkdir mask


def usage() :
   print("""

   To run basecalling:

      runAlbacore.py [-abD] [-j n] [-k kit] [-f flowcell_type]

      Must be called inside a run directory. Fast5 files will be searched
      recursively. Input and output are partitioned into n segments in the
      scratch directory.

         -j   number of jobs to use (default 40)
         -D   use 1D^2 chemistry
         -b   use the barcode keyword
         -a   get flowcell and kit from fast5 file
         -f   flowcell type (must specify kit too)
         -k   kit used in sequencing (must specify flowcell too)
      
      Output fastq is located under the stcratch/out directory. A RunInfo json
      is also created with all the experiment's details.

   To collect results:

      runAlbacore.py [-z ] -P
         
         -z   compress merged fastq

      Will merge fastq and summary files into one file. It also produces
      MinionQC for each summary file.
""")

def doBaseCalling(useFast5Info,flowcellType,kit,barcode,D2,chunks) :
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

   #initialize values
   flowcell = None
   start = None
   pc = None
   runNumber = None
   etype = None
   minionId = None
   minknowV = None

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
   deviceType =  f['/%s/tracking_id/'%globalKey].attrs['device_type'].decode()
   deviceId =  f['/%s/tracking_id/'%globalKey].attrs['device_id'].decode()
   minknowV =  f['/%s/tracking_id/'%globalKey].attrs['version'].decode()

   if ( D2 ) :
      chemistry = '1D^2'
   else :
      chemistry = '1D'

   #store all the run parameters in a json
   runInfo = { 'flowcell' : flowcell, 'flowcellType' : flowcellType_fast5 , 'kit': kit_fast5 , 'startTime' : start, 'albacoreVersion' : albacoreVersion , 'pc' : pc , 'runNumber' : runNumber, 'barcoded' : barcode, 'experimentType' : etype, 'instrumentType' : deviceType, 'instrumentId' : deviceId, 'minknowCoreVersion' : minknowV, 'chemistry' : chemistry }

   f = open('./RunInfo','w')
   json.dump(runInfo,f)
   f.close()

   if ( flowcell == None ) :
      print("Flowcell name is not defined in fast5, using None !")
   if (  flowcellType_fast5 == '' or kit_fast5 == '') :
     exit("Couldn't read run information from fast5 file (flowcell type %s, kit %s)!"%(flowcellType_fast5, kit_fast5) )

   if ( useFast5Info ) :
      flowcellType= flowcellType_fast5
      kit = kit_fast5
   else :
      if ( flowcellType != flowcellType_fast5 or kit != kit_fast5 ) :
         print("Warning specified flowcell %s and kit %s dont match fast5 info ( %s %s"%(flowcellType,kit,flowcellType_fast5,kit_fast5) )

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

#merge the resulting fastq and summary files from Albacore basecalling
#make MinionQC plots for each summary file
#compress the merged fastqs
def doPost(pigz) :
   #use this pigz
   pigz = '/apps/PIGZ/2.3.1/bin/pigz'

   #barcode names in the LIMS
   ONT2LIMS = {
            'barcode01' : 'NB01',
            'barcode02' : 'NB02',
            'barcode03' : 'NB03',
            'barcode04' : 'NB04',
            'barcode05' : 'NB05',
            'barcode06' : 'NB06',
            'barcode07' : 'NB07',
            'barcode08' : 'NB08',
            'barcode09' : 'NB09',
            'barcode10' : 'NB10',
            'barcode11' : 'NB11',
            'barcode12' : 'NB12',
            'barcode13' : 'NB13',
            'barcode14' : 'NB14',
            'barcode15' : 'NB15',
            'barcode16' : 'NB16',
            'barcode17' : 'NB17',
            'barcode18' : 'NB18',
            'barcode19' : 'NB19',
            'barcode20' : 'NB20',
            'barcode21' : 'NB21',
            'barcode22' : 'NB22',
            'barcode23' : 'NB23',
            'barcode24' : 'NB24',
            'unclassified' : 'unclassified'
            }

   #look for the RunInfo file generated by the basecaller script
   if ( not os.path.isfile('./RunInfo') ) :
      sys.exit( 'File RunInfo not found! Must run basecalling script to generate it.\n')
   f= open('./RunInfo','r')
   runInfo = json.load(f)
   f.close()

   #check the barcodes
   if ( runInfo['barcoded'] ) :
      cmd =  "lfs find ./scratch -type d -name 'barcode*'  2> /dev/null "
      p = Popen(cmd ,shell=True, stdout=PIPE)
      (output, error) = p.communicate()
      p.wait()
      fout = output.decode().split('\n')[:-1]         
      barcodes = set()
      for iout in fout :
         ibar = iout.split('/')[-1]
         if ( ibar not in ONT2LIMS ) : continue #barcode must be in the LIMS otherwise ignore
         barcodes.add( ibar )

   #setup script 
   outfile = open('.' + '/post.cmd','w')
   outfile.write(
   """#!/bin/bash
   # @ output = post.out
   # @ error = post.err
   # @ total_tasks = 1
   # @ cpus_per_task = 4
   # @ wall_clock_limit = 02:59:00
   """)
   outfile.write("# @ initialdir = .\n")
   outfile.write("# @ job_name = %s_%s_post\n"%(runInfo['flowcell'],1 ))
   outfile.write("\n")
   outfile.write("module purge\n")
   outfile.write("module load gcc/6.3.0\n")
   outfile.write("module load R/3.5.0\n")
   outfile.write("export R_LIBS_USER=\"/home/sbsuser/MinION/R\"\n\n")
   outfile.write("set -e\n\n")

   #merge all the fastq files (exclude 1D2 reads)
   if (  runInfo['barcoded'] ) :
      for ibarcode in ( list(barcodes) + ['unclassified'] )  :
         fliName = "%s_%s_%s"%(runInfo['flowcell'],runInfo['runNumber'], ONT2LIMS[ibarcode])
         fastqFile = "%s/%s_%s_%s.fastq"%('.',runInfo['flowcell'],runInfo['runNumber'], ONT2LIMS[ibarcode] )            
         outfile.write("lfs find %s -name '*.fastq'  | grep -v '/1dsq_analysis/' | grep '/%s/' | xargs -i cat {} > %s\n"%('.'+"/scratch/out",ibarcode,fastqFile) )
   else :
      fliName = "%s_%s_%s"%(runInfo['flowcell'],runInfo['runNumber'], '0')
      fastqFile = "%s/%s_%s_0.fastq"%('.',runInfo['flowcell'],runInfo['runNumber'])
      outfile.write("lfs find %s -name '*.fastq' |  grep -v '/1dsq_analysis/' | xargs -i cat {} > %s\n"%('.'+"/scratch/out",fastqFile) )

   #merge all the fastq files for 1D2 reads
   if ( runInfo['chemistry'] == '1D^2' ) :
      if ( not os.path.exists('./1dsq_analysis') ) : os.makedirs('./1dsq_analysis')
      fastqFile = "%s/%s_%s_0.fastq"%('./1dsq_analysis',runInfo['flowcell'],runInfo['runNumber'])
      outfile.write("lfs find %s -name '*.fastq' |  grep '/1dsq_analysis/' | xargs -i cat {} > %s\n"%('.'+"/scratch/out",fastqFile) )

   #merge all the sequencing summary files (exclude 1D2 reads) 
   if (  runInfo['barcoded'] ) :
      for ibarcode in ( list(barcodes) + ['unclassified'] )  :
         outfile.write("lfs find %s/scratch/out -name 'sequencing_summary.txt' | grep -v '/1dsq_analysis/'  | head -1 | xargs head -1 > %s/sequencing_summary.%s.txt\n"%('.','.', ONT2LIMS[ibarcode]))
         outfile.write("lfs find %s/scratch/out -name 'sequencing_summary.txt' | grep -v '/1dsq_analysis/'  | xargs -i awk \'$20 == \"%s\"\' {}  >> %s/sequencing_summary.%s.txt\n"%('.',ibarcode,'.', ONT2LIMS[ibarcode]))
         #generate MinionQC plots
         outfile.write("Rscript /home/sbsuser/MinION/minion_qc-1.3.0/MinionQC.R -f %s -i %s -o %s/qc_plots/%s\n"%(fliName,'.'+'/sequencing_summary.%s.txt'%ONT2LIMS[ibarcode],'.', ONT2LIMS[ibarcode] ))
   else :
      outfile.write("lfs find %s/scratch/out -name 'sequencing_summary.txt' | grep -v '/1dsq_analysis/' | head -1 | xargs head -1 > %s/sequencing_summary.0.txt\n"%('.','.'))
      outfile.write("lfs find %s/scratch/out -name 'sequencing_summary.txt' | grep -v '/1dsq_analysis/' | xargs -i tail -n+2 {}  >> %s/sequencing_summary.0.txt\n"%('.','.'))
      outfile.write("Rscript /home/sbsuser/MinION/minion_qc-1.3.0/MinionQC.R -f %s -i %s -o %s/qc_plots/0\n"%(fliName,'.'+'/sequencing_summary.0.txt','.'))

   #merge all the sequencing summary files for 1D2 reads
   if ( runInfo['chemistry'] == '1D^2' ) :
      outfile.write("lfs find %s/scratch/out -name 'sequencing*summary.txt' | grep '/1dsq_analysis/' | head -1 | xargs head -1 > %s/sequencing_summary.txt\n"%('.','1dsq_analysis'))
      outfile.write("lfs find %s/scratch/out -name 'sequencing*summary.txt' | grep '/1dsq_analysis/' | xargs -i tail -n+2 {}  >> %s/sequencing_summary.txt\n"%('.','1dsq_analysis'))

   #compress all the merged fastq files
   if ( zz ) :
      outfile.write("%s -p4 *.fastq\n"%pigz)
      if ( runInfo['chemistry'] == '1D^2' ) :
         outfile.write("%s -p4 ./1dsq_analysis/*.fastq\n"%pigz)

   #touch post done to show we completed successfully 
   outfile.write("touch post.done\n\n")
   outfile.close()
   launch = ['/opt/perf/bin/mnsubmit', '.' + '/post.cmd' ]
   p = Popen(launch, stdout=PIPE, stderr=PIPE)
   p.wait()
   out = p.stdout.read().decode()
   err = p.stderr.read().decode()
   if ('ERROR' in err) :
      print("!! Error submitting job array to cluster (%s)"%(err))
      raise Exception("Error submitting job array to cluster")
   print(out)




try:
   opts, args = getopt.getopt(sys.argv[1:],"k:f:j:abDPz", ["help"])
except getopt.GetoptError as err:
# print help information and exit
   print(str(err))
   usage()
   exit()

kit = None
flowcellType = None
barcode = False
useFast5Info = False
D2 = False
chunks = 40
post = False
zz = False

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
   if iopt == "-P" :
      post = True
   if iopt == "-z" :
      zz = True


if ( post ) :
   doPost(zz)
   sys.exit(0)

if (  kit == None and flowcellType == None and useFast5Info == False) :
   usage()
   exit()

if ( ( kit != None or flowcellType != None) and useFast5Info == True ) :
   print("cant specify option -a and manually select flowcell and kit")
   usage()
   exit()

if ( ( kit == None or flowcellType == None) and useFast5Info == False ) :
   print("must specify both flowcell and kit")
   usage()
   exit()

doBaseCalling(useFast5Info,flowcellType,kit,barcode,D2,chunks)
