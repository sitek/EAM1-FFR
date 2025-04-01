
"""
This script processes EEG data to find trigger lags.

The script analyzes stimulus channels and correlates them
with stimulus audio files.
The results are saved to a CSV file.
"""

import os
import mne
import pandas as pd
from scipy.signal import butter, filtfilt, hilbert
from scipy.io import wavfile
import numpy as np

# Directories and paths
maindir = os.path.join('/Users/dsj3886/Library/CloudStorage/',
                       'OneDrive-SharedLibraries-NorthwesternUniversity/',
                       'SoundBrain Lab - Documents')
expdir = os.path.join(maindir, 'Lab Research Projects',
                      'Experiments/NU_Experiments/EAM1')
bidsdir = os.path.join(expdir, 'data/EEG/data-bids')
stimdir = os.path.join(expdir, 'K01_FFR/button_FFR/')
outname = os.path.join(expdir, 'data/EEG/EHL1_adjusted_triggertimes.csv')

# Read completed data
comp = pd.read_csv(outname)
fdone = comp['filename'].unique()

trigtimes = []

# Read subject folders
subfolders = [f for f in os.listdir(bidsdir) if f.startswith('sub')]

for subfolder in subfolders:
    subpath = os.path.join(bidsdir, subfolder, 'eeg')
    subfiles = [f for f in os.listdir(subpath) if f.endswith('.bdf')]

    for fn in subfiles:
        if fn in fdone:
            continue

        # Placeholder for reading EEG data (replace with appropriate function)
        # EEG = read_eeg_file(os.path.join(subpath, fn))

        # Extract stimulus channel
        # stimchan = EEG['data']['Erg1']

        # MNE version
        EEG = mne.io.read_raw_bdf(fn, preload=True)
        stimchan = EEG.pick(['Erg1'])[0]

        # Extract triggers and remove 255
        timestamps = np.array([event['latency'] for event in EEG['event']
                               if event['type'] != 255])

        # Filter stimulus channel
        b, a = butter(1, 1 * 2 / EEG['srate'], 'high')
        stim = filtfilt(b, a, stimchan)
        stim /= np.max(np.abs(stim))

        # Get stimulus wav files
        namepts = fn.split('_')
        order = int(namepts[3][4])
        task = namepts[2]

        if task in ['task-alice', 'task-mix']:
            stimwavs = (list(range(1, 16))
                        if order == 1 else list(range(16, 31)))
        elif task == 'task-einstein':
            stimwavs = list(range(1, 16))

        for tracknow, tracknum in enumerate(stimwavs, 1):
            if task in ['task-alice', 'task-mix']:
                wav_path = os.path.join(stimdir, 
                                        f'alice_stimuli/track{tracknum}.wav')
            else:
                wav_path = os.path.join(stimdir, 
                                        f'einst_stimuli/Track{tracknum}.wav')

            fs, y = wavfile.read(wav_path)

            # Get envelope of wav file
            stimsound = np.abs(hilbert(y))
            stimsound = np.interp(np.arange(0, 
                                            len(stimsound), 
                                            fs / EEG['srate']),
                                  np.arange(len(stimsound)), 
                                  stimsound)
            stimsound /= np.max(np.abs(stimsound))

            # Calculate delay
            corr = np.correlate(stimsound, stim, 'full')
            del_samples = np.argmax(corr) - (len(stimsound) - 1)

            del2 = (del_samples - timestamps[tracknow - 1]) / EEG['srate']
            samps = round(del2 * 128)

            trigtimes.append([fn, tracknow, del2, samps])

# Save adjusted trigtimes
headers = ['filename', 'trigger', 'newtime_sec', 'samps_to_add']
Tt = pd.DataFrame(trigtimes, columns=headers)

# Append new trig times to old csv file
T = pd.concat([comp, Tt])
T.to_csv(outname, index=False)
