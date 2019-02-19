# ONT utilities
Various utilities for analyzing Oxford Nanopore experiments

### runAlbacore.py
Basecall an ONT experiment using the Albacore basecaller (no longer being mantained). The script must be invoked in the directory tree that contains the fast5 to be processed, all subdirectories are searched recursively for input.

A scratch directory is created with multiple input and output subfolders (one for each job). A Slurm job array is created and submitted to the execution queue with output sequencing_summary and fastq files located in the out directories. Each Albacore instance is asigned 4 threads (anything beyond that doesn't seem to have much of an impact).

Demultiplexing is also carried out at this step for barcoded experiments.

A RunInfo json file is generated with all the details of the experiment and basecalling.

    runAlbacore.py [-abD] [-j n] [-k kit] [-c conf]  [-f flowcell_type]

      Create chunks of raw fast5 input data to run several instances of 
      the Guppy basecaller in a cluster.

         -j   number of jobs to use (default 40)
         -D   use 1D^2 chemistry
         -b   use the barcode keyword
         -a   get flowcell and kit from fast5 file
         -f   flowcell type (must specify kit too)
         -k   kit used in sequencing (must specify flowcell too)

      Input fast5 must be located in the current directory tree.
      
      Output fastq is located under the stcratch/out directory. A RunInfo json
      is also created with all the experiment's details
      

### runGuppy.py
Basecall an ONT experiment using the Guppy basecaller. The script must be invoked in the directory tree that contains the fast5 to be processed, all subdirectories are searched recursively for input.

A scratch directory is created with multiple input and output subfolders (one for each job). A Slurm job array is created and submitted to the execution queue with output sequencing_summary and fastq files located in the out directories. Each Guppy instance is assigned 8 threads ( one devoted reading thread ).

If basecalling a multiplexed experiment an additional demultiplexing step must be done after this.

A RunInfo json file is generated with all the details of the experiment and basecalling.

    runGuppy.py [-aD] [-j n] [-k kit] [-c conf]  [-f flowcell_type]

      Create chunks of raw fast5 input data to run several instances of 
      the Guppy basecaller in a cluster.

         -j   number of jobs to use (default 30)
         -D   use 1D^2 chemistry
         -a   get flowcell and kit from fast5 file
         -f   flowcell type (must specify kit too)
         -k   kit used in sequencing (must specify flowcell too)
         -c   specify a configuration file (instead of flowcell & kit)

      Input fast5 must be located in the current directory tree.
      
      Output fastq is located under the stcratch/out directory. A RunInfo json
      is also created with all the experiment's details


### runPost.py
Final step in the bascalling process. The script merges fastq and sequencing_summary files. It also produces a very nice set of plots to asses quality of a run using the MinIONQC R package. It also compreses the resulting fastq files using the pigz command.

The script submits a short 4 thread job to the execution queue trough Slurm.

Script must be called inside of the run directory. Runs must have been demultiplexed prior to calling this script.
