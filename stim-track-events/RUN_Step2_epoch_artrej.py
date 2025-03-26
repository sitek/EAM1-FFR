#This script was written by Jacie R. McHaney for analyzing continuous speech EEG data.
#this is STEP 2 in the process.
#This script will epoch the data and run artifact rejection
#created on 12/2/2024

import os
import mne
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import pandas as pd


# directories
main_dir = os.path.join('/projects/b1208/EHL1')
mix_dir = os.path.join(main_dir, 'data/EEG')
data_dir = os.path.join(mix_dir, 'processed/ref_down_filt')
out_dir = os.path.join(data_dir, 'epoch_artrej')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Parameters
newsr = 128
trigs = list(range(1,4))
epochmin = -2
epochmax = 70 #shortest alice track in seconds
montname = os.path.join(main_dir, 'analysis/EEG/biosemi_32ch_2mastoid_locs.csv')
trigtimes = pd.read_csv(f"{main_dir}/data/EEG/EHL1_adjusted_triggertimes.csv") #get this from `FIRSTSTEP_findtriglags.m` script in matlab.
n_evts = 15

#get list of files processed in step 1
files = [i for i in os.listdir(data_dir) if i.endswith('.fif')]

#go file by file and do the following:
#1. epoch the data
#2. run artifact rejection on the epoched data

for f in range(len(files)):
    
    #read in current subject file
    fn = files[f]
    fif_path = os.path.join(data_dir, fn)
    data_filtered = mne.io.read_raw_fif(fif_path, preload=True)

    #find trigger events
    events = mne.find_events(data_filtered, stim_channel = 'Status', initial_event=False)

    #remove unwanted trigger event 255
    events = mne.pick_events(events, exclude=255)
    print(events)

    #modify trigger times
    namepts = fn.split('_')
    newfn = f'{namepts[0]}_{namepts[1]}_{namepts[2]}_{namepts[3]}'
    newevts = events
    subtimes = trigtimes[trigtimes['filename'].str.contains(newfn)]

    for evt in range(0,n_evts):
        tmptrig = evt+1
        currinfo = subtimes[subtimes['trigger']==tmptrig]
        trigadj = currinfo['samps_to_add'].tolist()
        newevts[evt,0]=newevts[evt,0]+trigadj

    ### Epoch data
    epochs = mne.Epochs(data_filtered, 
                    events = newevts, 
                    tmin = epochmin, tmax = epochmax, 
                    baseline = None,
                    preload=True)
    
    ### Artifact Rejection
    # Perform regression using the EOG sensor as independent variable and the EEG sensors as dependent variables.
    model_plain = mne.preprocessing.EOGRegression(picks="eeg", picks_artifact="eog").fit(epochs)

    #now substract EOG from the EEG
    epochs_clean_plain = model_plain.apply(epochs)

    # After regression, we should redo the baseline correction
    epochs_clean_plain.apply_baseline()

    # create epochs with the evoked subtracted out
    epochs_sub = epochs.copy().subtract_evoked()

    # perform regression
    model_sub = mne.preprocessing.EOGRegression(picks="eeg", picks_artifact="eog").fit(epochs_sub)

    # apply the regression coefficients to the original epochs
    epochs_clean_sub = model_plain.apply(epochs).apply_baseline()

    ### save
    cleandat = epochs_clean_sub
    Toutname = fn[0:len(fn)-4]
    outname = Toutname + "_epoch_artrej.fif"
    savename = os.path.join(out_dir, outname)
    cleandat.save(savename)