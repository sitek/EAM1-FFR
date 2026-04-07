# Auditory–motor interactions in the early auditory system
Frequency-following response (FFR) analysis code for the EAM1 project

### Manuscript details
Manuscript in preparation – stay tuned for a preprint, coming soon!

### Data availability
**TODO**: Share data on [OpenNeuro](https://openneuro.org)

### Analysis steps
1. Extract precise sound onset timing from the EEG ergo channel: [stim-track-events/preprocess_data.m](stim-track-events/preprocess_data.m)
  - **TODO**: Update the `stim-track-events` code to a modern Python implementation (in progress)
2. Convert the EEG data to the Brain Imaging Data Structure (EEG-BIDS): [bids-conversion.ipynb](bids-conversion.ipynb)
3. Preprocess and analyze the FFRs: [ffr_processing-minimal_mne-bids.ipynb](ffr_processing-minimal_mne-bids.ipynb)
4. Preprocess and analyze the cortical ERPs: [erp_processing-minimal_mne-bids.ipynb](erp_processing-minimal_mne-bids.ipynb)


## How to contribute
Windows: in a bash terminal, run: `bash setup-git-jupyter-filters.sh`
Mac: run `chmod +x setup-git-filters.sh` then `./setup-git-filters.sh `