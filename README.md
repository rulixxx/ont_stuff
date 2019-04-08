# ONT utilities
Various utilities for efficiently basecalling Oxford Nanopore experiments in an HPC environment.

### runAlbacore.py
Basecall an ONT experiment using the Albacore basecaller (no longer being mantained). The script must be invoked in the directory tree that contains the fast5 to be processed, all subdirectories are searched recursively for input.

A scratch directory is created with multiple input and output subfolders (one for each job). A Slurm job array is created and submitted to the execution queue with output sequencing_summary and fastq files located in the out directories. Each Albacore instance is assigned 4 threads (anything beyond that doesn't seem to have much of an impact).

The complete processing of an experiment requires two steps: 1) basecalling & 2) collection of results.

Demultiplexing is carried out for barcoded experiments if the barcoding flag is specified.


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

Scratch directory should be removed upon completion to reduce inode usage.
```
lfs find ./scratch -type f | xargs rm
rm -rf ./scratch
```
   

### runGuppy.py
Basecall an ONT experiment using the Guppy basecaller ( CPU version). The script must be invoked in the directory tree that contains the fast5 to be processed, all subdirectories are searched recursively for input.

A scratch directory is created with multiple input and output subfolders (one for each job). A Slurm job array is created and submitted to the execution queue with output sequencing_summary and fastq files located in the out directories. Each Guppy instance is assigned 8 threads ( one devoted reading thread ).

If basecalling a multiplexed experiment an additional demultiplexing step must be done after basecalling. Demultiplexing is carried out using qcat ( https://github.com/nanoporetech/qcat ). After comparison with other demultiplexers it was deemed the best performer. (https://community.nanoporetech.com/posts/demultiplexer-comparison)

After basecalling [demultiplexing] a post processing step generates merged fastq and summary files. It also runs the MinionQC ( https://github.com/roblanf/minion_qc ) R module to produce quality control plots.


    1) Run basecalling:   

      runGuppy.py [-aqD] [-j n] [-k kit] [-c conf]  [-f flowcell_type]

      Must be called inside a run directory. Fast5 files will be searched
      recursively. Input and output are partitioned into n segments in the
      scratch directory.

         -j   number of jobs to use (default 30)
         -D   use 1D^2 chemistry
         -a   get flowcell and kit from fast5 file
         -f   flowcell type (must specify kit too)
         -k   kit used in sequencing (must specify flowcell too)
         -c   specify a configuration file (instead of flowcell & kit)
         -q   apply quality filtering

      Output fastq is located under the stcratch/out directory. A RunInfo json
      is also created with all the experiment's details.

    2) Demultiplex fasts (using qcat) :
   
      runGuppy.py [-j n ] [-d directory ] -B

         -j   number of jobs to use (default 30)
         -d   directory where to look for fastqs ( default ./scratch/out )

    3) Collect results:

      runGuppy.py [-z ] -P
         
         -z   compress merged fastq

      Will merge fastq and summary files into one file. It also produces
      MinionQC for each summary file. 

Scratch directory should be removed upon completion to reduce inode usage.
```
lfs find ./scratch -type f | xargs rm
rm -rf ./scratch
```

### MinIONPip.py
Example of an automated basecalling pipeline for ONT experiments. It is in large part derived from the previous script but has been rewritten in order to improve fault tolerance, to allow complete automation and an integration with a LIMS.

The script carries out the following:

    * Check experiment with LIMS details
    * Decompress raw data
    * Split input reads into batches to allow parallel processing
    * Submit a job array to do the basecalling with Guppy
    * Demultiplex fastqs with qcat using a job array
    * Compute stats and post to the LIMS
    * Produce QC plots with MinionQC.R
    * Compress fastq files and remove temporary files

The pipeline is currently being run from a cron job and will look in the directories dataPC1 and dataPC2 for new experiments. It will scan for zip files and will start processing onece the file size remains constant. 

A nice featuree of the pipeline is that it reports most errors by email so that problems can be diagnosed and resolved quickly.

