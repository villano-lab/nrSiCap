## Current Release (v2.0.0)

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
