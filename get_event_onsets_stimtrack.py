"""
Processes raw EEG data (bdf format) from the experiment using MNE-Python and finds the onsets of the audio stimuli
(.wav) in the Erg1 channel.
Generates events.tsv files with stimulus onsets (seconds) and types (1 vs 2) for each participant, run, and condition.

To run script:
- Change 'top_dir' to the location of your data folder where the 'eeg_raw' and 'EAM2_stimuli' folders are located. (EAM2)
        or where data-bids folder is located (EAM1).
- Change 'save_dir' to the location where you want to save the generated events.tsv files. The last folder level must be
        'python' for the validation script to work.
- Change 'expected_peaks_num' to the number of stimuli you expect to see in each run. Default is 1200 for EAM1.
- Change 'skip_subjects' to an empty list if you want to run the script on all subjects. Otherwise, add the subject IDs
        you want to skip to the list.

Translated from the original MATLAB code by Laura 3/19/2026 (lauraraiff2030@u.northwestern.edu)
"""
import os
import mne
from scipy.signal import find_peaks, resample, correlate, correlation_lags, hilbert
from scipy.io import wavfile
import pandas as pd
import numpy as np


def normalized_xcorr(x, y):
    """ Calculates the normalized cross correlation between two signals."""
    corr = correlate(x, y, mode='full')
    norm = np.linalg.norm(x) * np.linalg.norm(y)

    if norm > 0:
        corr /= norm

    return corr

def finddelay(x, y, maxlag=None):
    """ 
    based on the matlab function finddelay: 
    
    The finddelay function uses the xcorr function to determine the cross-correlation between each pair of signals at all possible lags 
    specified by the user. The normalized cross-correlation between each pair of signals is then calculated. The estimated delay is given 
    by the negative of the lag for which the normalized cross-correlation has the largest absolute value. If more than one lag leads to the 
    largest absolute value of the cross-correlation, such as in the case of periodic signals, the delay is chosen as the negative of the 
    smallest (in absolute value) of such lags. 
    """

    # set maxlag value if not specified
    # if maxlag is None:
    #     [x_rows, x_cols] = np.shape(x)
    #     maxlag = max(x_rows, x_cols) - 1

    # first normalized cross correlation
    corr = normalized_xcorr(x, y)
    
    # estimated delay is negative of the lag for which the normalized cross correlation has the largest absolute value
    lags = correlation_lags(len(x), len(y), mode='full')
    abs_corr = np.abs(corr)
    max_val = np.max(abs_corr)

    # handle ties: pick the smallest absolute lag, prefer positive/zero lag
    tied_lags = lags[abs_corr == max_val]
    pos_tied = tied_lags[tied_lags >= 0]
    neg_tied = tied_lags[tied_lags < 0]

    if len(pos_tied) > 0 and len(neg_tied) > 0:
        best_pos = pos_tied[np.argmin(pos_tied)]  # smallest positive
        best_neg = neg_tied[np.argmin(np.abs(neg_tied))]  # smallest abs negative
        best_lag = best_pos if best_pos < abs(best_neg) else best_neg
    elif len(pos_tied) > 0:
        best_lag = pos_tied[np.argmin(pos_tied)]
    else:
        best_lag = neg_tied[np.argmin(np.abs(neg_tied))]

    return -best_lag


def envelope(x):
    # remove DC offset
    x_mean = np.mean(x)
    x_centered = x - x_mean

    # envelope amplitude
    x_amp = np.abs(hilbert(x_centered))

    # add offset back in 
    y_upper = x_mean + x_amp
    y_lower = x_mean - x_amp

    return y_upper, y_lower

# define file paths (CHANGE THESE TO YOUR FILE PATHS)
root_dir = r"C:\Users\Laura\OneDrive - Northwestern University\SoundBrain Lab - EAM1\data-bids"
save_dir = r"C:\Users\Laura\Documents\PhD\Soundbrain lab\EAM\data-bids\python"
expected_peaks_num = 1200
skip_subjects = ['sub-18']

all_sub_dirs = [d for d in os.listdir(root_dir) if d.startswith('sub-')]
for sub_id in all_sub_dirs:
    if sub_id in skip_subjects:
        continue
    # data paths
    stimulus_path = os.path.join(root_dir, 'stimuli')
    save_dir_sub = os.path.join(save_dir, sub_id)
    eeg_path = os.path.join(root_dir, sub_id, 'eeg')

    eeg_files = [f for f in os.listdir(eeg_path) if f.endswith('.bdf')]
    stimuli_files = [f for f in os.listdir(stimulus_path) if f.endswith('.wav')]

    # get eeg sampling frequency 
    try:
        first_eeg = mne.io.read_raw_bdf(os.path.join(eeg_path, eeg_files[0]), preload=False)
    except IndexError:
        print(f"No EEG files found for {sub_id} in {eeg_path}. Skipping this subject.")
        continue
    except Exception as e:
        print(f"Error occurred while reading EEG file: {e}")
        continue
    fs = first_eeg.info['sfreq']

    # load and resample stimulus audio files
    stimuli_audio = {}
    for file in stimuli_files:
        stimulus_file = os.path.join(stimulus_path, file)
        stim_fs, stim_file = wavfile.read(stimulus_file) 
        n_new_samples = int(len(stim_file) * fs / stim_fs)
        stimuli_audio[file] = resample(stim_file, n_new_samples)

    # init stim variables
    num_stimuli_arr = np.arange(1, len(stimuli_files) + 1)

    # load eeg data
    for file in eeg_files:
        # get information from file name
        eeg_file = os.path.join(eeg_path, file)
        fname_no_ext = file.split('.')[0]
        fname_split = fname_no_ext.split('_')
        sub = fname_split[0]
        cond_id = fname_split[1]
        run_id = fname_split[2]

        eeg_raw = mne.io.read_raw_bdf(eeg_file, preload=True)

        # find peaks in the stimtrack data
        stimtrack_data = eeg_raw.get_data(picks=['Erg1']).T.flatten()
        stimtrack_data = stimtrack_data / max(abs(stimtrack_data))  # normalize 

        # compute and normalize the upper envelope of the signal
        upper_env, _ = envelope(stimtrack_data)
        upper_env = upper_env / max(abs(upper_env))

        peaks, _ = find_peaks(upper_env, distance=fs*0.17, height=0.5)

        # we expect there to be 1200 peaks, but sometimes there are more or less
        if len(peaks) != expected_peaks_num:
            print(f"Warning: Expected 1200 peaks, but found {len(peaks)} for {file} with condition {cond_id} and run {run_id}")
            isi = np.diff(peaks)

            # for extra peaks in the beginning of the stimtrack
            if len(peaks) > expected_peaks_num:
                gap_idx = np.where(isi > np.median(isi) * 10)
                if gap_idx:
                    peaks = peaks[(gap_idx[0][0] + 1):]

        # epoching
        onset_arr = []
        stimuli_final = []
        stim_nums = []
        for i, peak in enumerate(peaks):
            epoch_quarter = round(0.15*fs)

            # check length, skipping if exceeds stim length or if epoch starts before 0
            epoch_idx = np.arange(peak - epoch_quarter, peak + epoch_quarter + 1, dtype=int)
            epoch_idx_max = max(epoch_idx)
            length_stim = len(stimtrack_data)

            if length_stim < epoch_idx_max:
                continue

            if (peak - epoch_quarter) < 0:
                continue
            
            stim_epoch = stimtrack_data[epoch_idx]

            # match epoch to stim
            max_corr = []
            max_corr_idx = []
            all_r_stim = []
            for s in range(len(stimuli_files)):
                stim = stimuli_audio[stimuli_files[s]]

                # pad the stimulus
                if len(stim) < len(stim_epoch):
                    stim_padded = np.pad(stim, (0, len(stim_epoch) - len(stim)), mode='constant')
                else:
                    stim_padded = stim[:len(stim_epoch)]

                # cross-correlate the stimulus with the epoch and find the index of the maximum correlation
                r_stim = normalized_xcorr(stim_epoch, stim_padded)
                max_corr.append(np.max(r_stim))
                max_corr_idx.append(np.argmax(r_stim))
                all_r_stim.append(r_stim)
            
            # get the type of stimulus 
            best_match_idx = np.argmax(max_corr)
            signal = stimuli_audio[stimuli_files[best_match_idx]]
            stimuli_final.append(signal)
            stim_nums.append(num_stimuli_arr[best_match_idx])

            # get the onset of the stimulus
            alignment_delay = finddelay(signal, stim_epoch)
            new_onset = peak - alignment_delay

            onset_arr.append(new_onset/fs)
        
        print(f"Generated {len(onset_arr)} onsets for {file} with condition {cond_id} and run {run_id}")

        # write to events.tsv file
        events_file = os.path.join(save_dir_sub, f"{sub_id}_{cond_id}_{run_id}_events.tsv")
        if not os.path.exists(save_dir_sub):
            os.makedirs(save_dir_sub)
        
        event_data = pd.DataFrame(columns=['type', 'onset'])
        event_data['type'] = pd.Series(stim_nums)
        event_data['onset'] = pd.Series(onset_arr)

        event_data.to_csv(events_file, sep='\t', index=False)

        print(f"Saved events file for {file} with condition {cond_id} and run {run_id} at {events_file}")

    print('done with subject')
