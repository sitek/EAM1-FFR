import os
import pandas as pd
from glob import glob
from mne_bids import BIDSPath, find_matching_paths

# Directory containing stimtrack files
stimtrack_dir = os.path.join('/Users/dsj3886/data_local/',
                             'EAM1_local/data-bids/derivatives/',
                             'events-stimtrack')
stimtrack_files = [f for f in os.listdir(stimtrack_dir) if f.endswith('.tsv')]
print(f'Found {len(stimtrack_files)} stimtrack files')

# Directory for BIDS-compatible events files
# This is where the new events files will be saved
stimtrack_bids_dir = os.path.join('/Users/dsj3886/data_local/',
                                  'EAM1_local/data-bids-stimtrack/')

def replace_events_file(stimtrack_file):
    """
    Create a new events file from a stimtrack-based events file.
    """
    stimtrack_df = pd.read_csv(stimtrack_file, sep='\t')

    # Extract the subject ID from the filename
    subject_id = os.path.basename(stimtrack_file).split('_')[0].replace('sub-', '')
    task_id = os.path.basename(stimtrack_file).split('_')[1].replace('task-', '')
    run_id = os.path.basename(stimtrack_file).split('_')[2].replace('run-', '')

    # Construct the corresponding events file path
    events_fname = f'sub-{subject_id}_task-{task_id}_run-{run_id}_events.tsv'
    events_file = os.path.join(stimtrack_bids_dir,
                               f'sub-{subject_id}', 
                               'eeg',
                               events_fname)
    
    events_df = stimtrack_df.copy()
    events_df.rename(columns={'type': 'value'}, inplace=True)
    events_df['duration'] = 0.170  # Set duration to 170 ms
    events_df['trial_type'] = 'unknown'  # Default value
    events_df['trial_type'][events_df['value'] == 1] = 'positive'
    events_df['trial_type'][events_df['value'] == 2] = 'negative'

    events_df = events_df[['onset', 'duration', 'trial_type', 'value']]

    events_df.to_csv(events_file, sep='\t', index=False)

    return events_file, events_df


# Iterate over all subjects and tasks to create events files
for fpath in glob(stimtrack_dir + '/sub-*'):
    subject_id = os.path.basename(fpath).split('_')[0].replace('sub-', '')
    print(f'Processing subject: {subject_id}')
    for task_id in ['active', 'passive']:
        for run_id in ['1', '2', '3', '4', '5', '6']:
            events_fname = f'sub-{subject_id}_task-{task_id}_run-{run_id}_stimtrack_events.tsv'
            stimtrack_file = os.path.join(stimtrack_dir,
                                          events_fname)
            print(stimtrack_file)
            if os.path.exists(stimtrack_file):
                events_file, events_df = replace_events_file(stimtrack_file)
                print(f'Created events file: {events_file}')
