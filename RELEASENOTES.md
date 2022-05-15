## Next Release (v2.1.1)

## Release (v2.1.0) Date: 22.05.14

* Added a .ipynb notebook to calculate the final result for the paper only. It is in
  `1-PRD_Figures` and it's called `plot_final_yield.ipynb`. Note that it needs two large data
files from the OSF link: `mcmc_Sor_128walk_50kstep_SNorm_v3.pkl` and `mcmc_Sor_128walk_50kstep_SNorm_v4.pkl`

## Release (v2.0.0) Date: 21.12.14

Reorganized the repository:
* The `analysis_notebooks` directory has been removed. 
	* Any code kept from this directory has been moved to `0-Analysis/processing`.
* Python scripts with duplicates stored in multiple places are now only stored in the top-level `python` directory.
* Removed the directory `1-PRL_Figures` and replaced it with `1-PRD_Figures`.
	* Added notebooks for figures in PRD draft that were not in PRL draft.
	* Renumbered figure notebooks to match the new paper's order.
Documentation:
* Updated the README to match reorganization.
* Added instructions for retrieving large data files from OSF.io to the README.
* Added `CONTRIBUTING.md`.

Other Changes:
* Created a script `retrieve_data.sh` for easier retrieval of necessary data files from OSF.io.
* Pruned the environment in `nrSiCap.yml` to speed up environment creation time.
* Implemented testing through Travis-CI.



## Release (v1.0.1) Date: 21.10.27

* (hotfix) To update top name in README.md file

## Release (v1.0.0) Date: 21.10.23

* First release with very basic code working and going public. 
