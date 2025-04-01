%this script reads in BDF files and finds the lag between trigger and
%stimulus onset for ehl1
%created by Jacie R. McHaney on 12/6/24

close all
clear all

%directories and paths.
maindir = '/Users/dsj3886/Library/CloudStorage/OneDrive-SharedLibraries-NorthwesternUniversity/SoundBrain Lab - Documents';
addpath ([maindir '/Jacie/Pitt_CRC/eeglab14_1_2b']);
eeglab;

expdir  = [maindir '/Lab Research Projects/Experiments/NU_Experiments/EAM1'];
eegdir  = [expdir '/data-bids'];
stimdir = [expdir  'K01_FFR/button_FFR/stimuli/'];

%pull in already completed data
outname = [expdir 'data/EEG/EHL1_adjusted_triggertimes.csv'];
comp = readtable(outname);
fdone = unique(comp.filename); %get completed file names

%read in subject folders
subfolders = dir([rawdir '/sub*']);
subfolders = cellstr(char(subfolders.name));

trigtimes = {};

for i =1:length(subfolders)

    %get subject folders
    subpath = [rawdir '/' subfolders{i}];
    subfiles = dir([subpath '/eeg/*.bdf']);
    subfiles = cellstr(char(subfiles.name));

    for j = 1:length(subfiles)
        
        %read in EEG file
        EEG = [];
        fn = subfiles{j};

        %check if this file has already been adjusted and skip if so
        x = find(ismember(fdone,fn));
        if ~isempty(x)
            continue
        end

        %now read in data file
        EEG = pop_biosig([subpath '/' fn]);

        %extract stimulus channel
        chaninfo = struct2cell(EEG.chanlocs);
        chaninfo = cellstr(char(chaninfo(1,:)));
        idx = find(ismember(chaninfo,'Erg1'));
        stimchan = EEG.data(idx,:);

        %extract triggers and remove 255
        marks = EEG.event;
        ixM = find([marks.type]==255);

        %get original trigger times
        timestamps = vertcat(EEG.event.latency); %latency of each marker
        timestamps(ixM) = []; %remove 255 latencies

        %find onset of stimtrack track
        dtstim = detrend(stimchan,1);
        stim = double(dtstim);
        [b,a] = butter(1,1*2/EEG.srate, 'high');
        stim = filtfilt(b,a,stim);
        stim =stim/max(abs(stim));

        %read in corresponding wave file of the stimulus
        namepts = strsplit(fn,'_');
        order = str2num(namepts{4}(5));
        task = namepts{3};

        %get track list of stim files based on task and acq order
        if strcmpi(task,'task-alice') || strcmpi(task,'task-mix')
            if order==1
                stimwavs = arrayfun(@num2str,1:15,'UniformOutput',0);
            elseif order ==2
                stimwavs = arrayfun(@num2str,16:30,'UniformOutput',0);
            % elseif order ==3
            %     stimwavs = arrayfun(@num2str,31:45,'UniformOutput',0);
            % elseif order==4
            %     stimwavs = arrayfun(@num2str,46:60,'UniformOutput',0);
            end
        elseif strcmpi(task,'task-einstein')
            stimwavs = arrayfun(@num2str,1:15,'UniformOutput',0);
        end

        %get number of tracks for this file
        ntracks = length(stimwavs);

        %read in stimulus wav files
        tracknow = 1;
        while tracknow < ntracks+1

            %read in appropriate stim file
            if strcmpi(task,'task-alice') || strcmpi(task,'task-mix')
                [y,fs] = audioread([stimdir '/alice_stimuli/track' stimwavs{tracknow} '.wav']);
            elseif strcmpi(task,'task-einstein')
                [y,fs] = audioread([stimdir '/einst_stimuli/Track' stimwavs{tracknow} '.wav']);
            end

            %get envelope of wav file and resample it
            stimsound = abs(hilbert(y));
            stimsound = resample(stimsound,EEG.srate,fs);
            stimsound =stimsound/max(abs(stimsound));
            
            %calculate delay in samples
            del = finddelay(stimsound,stim); %in samples at biosemi fs

            %calculate in s the delay
            del2 = (del-timestamps(tracknow))/EEG.srate;

            %samples to add at 128 Hz sampling rate (my script downsamples
            %to 128 Hz later on. If you don't resample, then make samps = del;
            samps = round(del2*128);

            temparr = {fn,tracknow,del2,samps};

            %plot to see if triggers are appropriately aligned now.
            % figure;
            % plot(stimsound,'k');
            % hold on
            % plot(stim,'b');
            % hold off

            %add together
            trigtimes = [trigtimes;temparr];

            tracknow = tracknow+1;
        end %while


    end %files

end %folders

%save adjusted trigtimes
headers = {'filename','trigger','newtime_sec','samps_to_add'};
Tt = cell2table(trigtimes,'VariableNames',headers);

%append new trig times to old csv file
T = [comp;Tt];
writetable(T,outname);

