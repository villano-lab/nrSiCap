# `nrSiCap` 
Analysis code and documentation for the [first measurement of isolated nuclear recoils](https://arxiv.org/abs/2110.02751).  

## Overview
This repo contains some of the analysis code and a little of the documentation for the analysis of UMN R68, which used a SuperCDMS Si HV prototype to search for low energy (n,gamma) NRs from a PuBe neutron source.

Much of the other importatnt analysis and code is (unfortunately) scattered across a few places. The [Run 68 Summary on the UMN Wiki](https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68) contains lots of operational details and some pieces of analysis development, particularly in the [Analysis Task List](https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:tasklist). The [K100 enote page](http://www.hep.umn.edu/cdms/cdms_restricted/K100.html) hosts a number of more formal analysis notes, especially those relating to cut development.

This repo mostly contains code to compare simulations with measured data and perform fitting/parameter estimation.

## Important Files
The `python` directory contains some sets of useful functions, many are self explanatory. 
* `R68_efficiencies.py` provides the finalized cut and trigger efficiencies for the measured data. The file contains links to the relevant analysis desciptions from which each value was obtained. 
* `R68_load.py` packages the messy business of loading measured and simulated datasets.
* `R68_yield.py` has paramaterized yield models. It defines a genralized `Yield` class which is an object that gets passed to different routines in the analysis.
* `R68_spec_tools.py` contains a number of methods for generating, manipulating, and plotting energy spectra. This is where simulated hits from Geant4 or (n,gamma) sims have the yield applied and get combined into measured events. This section has undergone the most development and there are still vestigial functions with similar names to newer, more useful functions. This should be cleaned up eventually, but the comments make clear which are defunct.

`R68_MCMC_MPI.py`, `R68_MCMC_process.ipynb`, and `R68_MCMC_plots_v2.ipynb` are the main trio used to run and evaluate MCMC fits. They are descibed more in the next section.

### Necessary Binaries
In order to run some of the code in this repository, you will need to retrieve some binary data files that were too large to be stored in the repository directly. 
These can be found on the [associated OSF project](https://osf.io/g4enq/). Download the following files and place them directly in the top-level `data` directory:
* Any files matching `mcmc_*_128walk_50kstep_SNorm_v?.pkl`
* `normsi_fast_200k.pkl`
* `v3_400k.pkl`
And within the directory `data/byseries/07180924_1710`, add the file `07180924_1710_traces.root` from the OSF directory `07180924_1710`.

### MCMC
The main workflow is to perform an MCMC fit using the script `0-Analysis/processing/R68_MCMC_MPI.py`. This runs an MCMC fit using [emcee](https://emcee.readthedocs.io/en/stable/) on a computing cluser. This requires a few things to be set up first. The user should be operating in the conda environment defined by r68_env.yml. They will also need to have installed [schwimmbad](https://schwimmbad.readthedocs.io/en/latest/install.html) and a compiled version of MPI like [OpenMPI](https://www.open-mpi.org/).

There are a bunch of settings in this file that define things like what data to load, which yield model to test, and how many MCMC steps to run for. They are all controlled by changing values in the `mcmc_data` dictionary that takes up the first ~100 lines of `R68_MCMC_MPI`.

Once everything is set up, the script can be run from the command line with something like:

> mpiexec --hostfile /home/mastx027/MPI/mpi_hostfile -n 24 python R68_MCMC_MPI.py

The reason this is organized to be run on a cluster is that it can take a long time to execute the number of MCMC trials to get a good amount of stats. Depending on the settings used, this can take from a few minutes to several days to complete, so I recommend running the command in a `screen` session that one can leave running and check back in on at a later time. When the calculations are complete, the scrupt will end by printing the filename where the results are stored. This name is automaticaly generated to be somewhat informative but also to not overwrite existing data, and is something of the form `data/mcmc_Sor_16walk_5000step_pois_v1.pkl`.

The notebook `0-Analysis/processing/R68_MCMC_process.ipynb` is designed to open an MCMC results file and do some basic analysis of the data including calculations of sampler autocovariance, sample thinning, parameter distributions, and best-fit spectra and yield distributions. These calculations are stored in the same file alongside the MCMC results. The notebook `0-Analysis/processing/R68_MCMC_plots_v2.ipynb` can then be run to make a set of standard plots of these results.

### Numbered Directories
There are three numbered directories - `0-Analysis`, `1-PRL_Figures`, and `2-Additional_Figures`. 
`0-Analysis` contains notebooks discussing steps of the analysis that may be of interest but were not critical enough to include in the paper ([arXiv](https://arxiv.org/abs/2110.02751)).
`1-PRL_Figures` contains code for generating any figures in the paper.
`2-Additional_Figures` contains code for generating additional figures. Some of these notebooks are placeholders, and others may overlap with code in the first two. These do not contain any notable annotation, interpretation, or description.

# Contributors  
Nick Mast  
Anthony Villano  
Kitty Harris  
Ana Gevorgian  

*Last updated 9 Nov 2021, v2.0.0*
