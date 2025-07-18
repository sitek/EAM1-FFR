%% Function to process the data

function preprocess_data(bidsroot, sub_id, out_dir)

%maindir = ['/Users/dsj3886/Library/CloudStorage/'...
%           'OneDrive-SharedLibraries-NorthwesternUniversity/'...
%           'SoundBrain Lab - Documents'];
%addpath ([maindir '/Jacie/Pitt_CRC/eeglab14_1_2b']);
eeglab;

% define stimulus files
stimpath = [bidsroot  '/stimuli/' ];

stim_id  = 'da'; 
stimfiles_cst = dir([stimpath, '*', stim_id,'*.wav']);
stimfiles = cellstr(char(stimfiles_cst.name));

% define EEG files
filenames = dir([bidsroot, '/', sub_id, '/eeg/*.bdf']);

% Loop through each element in the struct array
filePaths = cell(length(filenames), 1);
for i = 1:length(filenames)
    filePaths{i} = [filenames(i).folder '/' filenames(i).name];
end

% Display the result
disp(filePaths);

%% Load all files

%eeglab('redraw');

for i = 1:length(filePaths)
    try
        fname = filePaths{i};
    
        [~, basename, ~] = fileparts(fname);
        
        EEG = pop_biosig(fname); %, [], []);
        EEG = eeg_checkset( EEG );
    
        if strcmp(fname,'sub2039')==1
            continue
        end
        
        fname_split = strsplit(basename, '_'); 
        sub = fname_split{1};
        cond_id  =  fname_split{2};
        run_id   = fname_split{3};
    
        %% find peaks in EEG stimtrack
    
        fs = EEG.srate;
        
        chaninfo = struct2cell(EEG.chanlocs);
        chaninfo = cellstr(char(chaninfo(1,:)));
        idx = find(ismember(chaninfo,'Erg1'));
        stim = EEG.data(idx,:);
        
        [upperstimtrackenv, ~] = envelope(stim);
        %envlow = 10*delta;
        %[b,a] = butter(3,envlow, 'low');
        %upperstimtrackenv = filtfilt(b,a,upperstimtrackenv);
        %figure;plot(upperstimtrackenv)
        
        upperstimtrackenv = upperstimtrackenv/max(abs(upperstimtrackenv));
        stim =stim/max(abs(stim));
        
        [~, locs, ~] = findpeaks(upperstimtrackenv, ...
                                     'MinPeakDistance', 0.170*fs, ...
                                     'MinPeakHeight',0.2 );
        
        %figure;
        %plot(upperstimtrackenv), hold on, 
        %plot(stim, 'k'), plot(locs,peaks,'rx')
        
        %% stimfiles
        stimsound = {};
        for i = 1:length(stimfiles)
            
            [stimsound_tmp,stimfs] = audioread([stimpath, filesep, stimfiles{i}]);
            stimsound_tmp = resample(stimsound_tmp,fs,stimfs);
            stimsound{i} = stimsound_tmp;
        end
        %%
        stimcodes = 1:length(stimsound);
        stimulus = [];
    
        ffrepochs = [];
        ffrepochs_resp = [];
        
        %% epoching
    
        onset_arr = [];
        for i = 1:length(locs)
            
            quarterepoch = round(0.15*fs);
            %quarterepoch_resp = round(0.175*resp_fs);
        
            % check length
            ln1 = locs(i)-quarterepoch:locs(i)+quarterepoch;
            ln1_max = max(ln1);
            lnstim = length(stim);
            if (lnstim < ln1_max) == 1
                continue
            end
            
            tempstim = stim(locs(i)-quarterepoch:locs(i)+quarterepoch);
    
            maxr = []; 
            maxloc =[]; 
            all_Rstim = [];
            for k = 1:length(stimsound)
                
                tempsound = stimsound{k};
                zerolen = zeros(length(tempstim) - length(tempsound),1);
                tempsound = [tempsound ; zerolen];
                [Rstim, ~] = xcorr(tempsound,tempstim, 'coeff');
                [maxr(k), maxloc(k)] = max(Rstim);
                all_Rstim(:,k) = Rstim;
            end
        
            [~, matchdstimloc] = max(maxr);
        
            stimulus(i) = stimcodes(matchdstimloc);
            
            signal = stimsound{matchdstimloc};
            onsetcut = finddelay(signal,tempstim);
            
            new_onset = locs(i)-onsetcut;
            onset_arr(i) = new_onset/fs;
        
        end 
        
        nepochs = length(onset_arr);
        disp(['generated ', num2str(nepochs), ' epochs'])
    
        
        %% TODO: save as events.tsv file
    
        outname = strjoin({sub, cond_id, run_id, ...
                          'stimtrack', 'events.tsv'}, '_');
        outpath = [out_dir '/' outname];
    
        % Define the column headers and data array
        events_hdrs = {'type', 'onset'};
    
        events_data = cell(1200, 2);
        for i = 1:1200
            events_data{i, 1} = [num2str(stimulus(i))];
            events_data{i, 2} = [num2str(onset_arr(i))];
        end
        
        % Open a file for writing
        fileID = fopen(outpath, 'w');
        
        % Write the column headers
        fprintf(fileID, '%s\t%s\n', events_hdrs{:});
        
        % Write the data array
        for row = 1:size(events_data, 1)
            fprintf(fileID, '%s\t%s\n', events_data{row, 1}, events_data{row, 2});
        end
        
        % Close the file
        fclose(fileID);
        
        disp(['Data has been written to ', outpath]);
    catch
        disp(['error processing ', fname])
        continue;
    %eeglab redraw
    
    
    end
end