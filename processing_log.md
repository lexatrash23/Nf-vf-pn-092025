# Processing Log

This document details what I have done to troubleshoot issues with the venomflow pipeline.

## TO DO

[ ] Fix interproscan and signalP
[ ] Don't use combine, use metadata and tuples with join
[ ] Change input to be able to take multiple samples 
[ ] Make the process names more descriptive

## 2025-12-02 - Lei

FASRC has blocked the default anaconda channels from the cluster due to their potentially stopping free-use of the default channels. (see https://docs.rc.fas.harvard.edu/kb/python-package-installation/) Removed the following channels from yaml files after timeout issues creating environment:

```
  - defaults
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
```

Error:
```
      LibMambaUnsatisfiableError: Encountered problems while solving:
        - nothing provides __archspec 1 zen2 needed by _x86_64-microarch-level-3-2_zen2
      
      Could not solve for environment specs
      The following package could not be installed
      └─ _x86_64-microarch-level ==3 2_zen2 is not installable because it requires
         └─ __archspec ==1 zen2, which is missing on the system.
```
Deleted lines from the busco.yaml file:

```
  - _x86_64-microarch-level=3=2_zen2
```

Sometimes the process will just hang. Why?

```
Dec-03 11:20:33.673 [Task monitor] DEBUG n.processor.TaskPollingMonitor - !! executor slurm > No more task to compute -- The following nodes are still active:
```

## 2025-12-03 - Lei

The conda environment continues to fail to be created. I'm going to try changing the dependencies for tasks that using containers if available.

All dependencies have been changed to use containers with some exceptions. 

I used the command `nextflow inspect -concretize main.nf -profile cannon` to pre-download and test the existence for all container images. 

**Note** There may be some non-optimal containers. I didn't exhaustively search for the "official version" of every software. So things from biocontainers might have better/more stable versions. 

**question**: Is `combine` the right operator to use in these channel manipulations? I think it probably isn't, because it creates the Cartesian product of the channels.

PostTrimFastqc killed 137 due to OOM. Added default values for memory and time to nextflow.config.

Need to edit the `signalP`, `interproscan` processes to use the bin executables. 

## 2025-12-05 - Lei

Added symlnks to the signalp and interproscan software to the bin. SignalP is working, which presumably means interproscan will also work. 

BLASTX is failing:

```
Work dir:
  /n/holylfs05/LABS/informatics/Everyone/support/20251202-NaiduPraveena-Venomflow/venomflow/work/f7/e5625892f0312f31efdb75389c3e60

Command error:
  BLAST Database error: No alias or index file found for protein database [unitox_curated] in search path [/tmp/nxf.4VffydRGHg:/blast/blastdb:/blast/blastdb_custom:]. Please verify the spelling of the BLAST database and its molecule type.
```

Also the busco container is failing to be mounted properly

```
WARNING: skipping mount of /lineages: stat /lineages: no such file or directory
FATAL:   container creation failed: mount /lineages->/lineages error: while mounting /lineages: mount source /lineages doesn't exist
cp: cannot stat '.command.trace': No such file or directory
```

