FFR analysis code for EAM1 project

Steps:
1. Extract precise sound onset timing from the EEG ergo channel: [stim-track-events/preprocess_data.m](stim-track-events/preprocess_data.m)
2. Convert the EEG data to the Brain Imaging Data Structure (EEG-BIDS): [bids-conversion.ipynb](bids-conversion.ipynb)
3. Preprocess and analyze the FFRs: [ffr_processing-minimal_mne-bids.ipynb](ffr_processing-minimal_mne-bids.ipynb)
4. Preprocess and analyze the cortical ERPs: [erp_processing-minimal_mne-bids.ipynb](erp_processing-minimal_mne-bids.ipynb)
