%% changed the time points that are being analyzed for object project
% adjusted from original preprocessing script (run_preprocess) - Nov 8,
% 2023 AK

%% set up and set paths 

function shortened_preprocess(subjectnr)

    fprintf('preprocessing sub-%02i\n',subjectnr);tic; % prints in command window "preprocessing subject-##"

% cosmo and fieldtrip locations
    addpath('~/CoSMoMVPA/mvpa/')
    addpath('~/fieldtrip-20231015');
    ft_defaults;

% files and paths

    datapath = '/Users/f004p7b/Documents/australia/objects_project/analysis_bids/data'; %where raw data is located
    mkdir(sprintf('%s/derivatives/cosmomvpa',datapath)); %makes derivatives folder (bids format)
    
    cosmofn = sprintf('%s/derivatives/cosmomvpa/sub-%02i_cosmomvpa.mat',datapath,subjectnr);
    behavfn = sprintf('%s/sub-%02i/eeg/sub-%02i_task-targets_events.tsv',datapath,subjectnr,subjectnr);
    source_behavfn = sprintf('%s/sourcedata/sub-%02i_task-targets.csv',datapath,subjectnr);
    source_filename = sprintf('%s/sourcedata/sub-%02i_task-targets_eeg.bdf',datapath,subjectnr);
    raw_filename = sprintf('%s/sub-%02i/eeg/sub-%02i_task-targets_eeg.bdf',datapath,subjectnr,subjectnr);
 
    if ~exist(raw_filename,'file')
        mkdir(fileparts(raw_filename))
        movefile(source_filename,raw_filename)
    end

    %% read events 

    hdr = ft_read_header(raw_filename); %loads the raw EEG data

    T = readtable(source_behavfn); % loads the behavioral info/data about which trials occurred when & the stimulus that was presented 

    % find the STATUS channel and read the values from it - I think the
    % status channel is the channel where signals are sent from experiment
    % to the EEG machine
    schan = find(strcmpi(hdr.label,'STATUS'));
    sdata = ft_read_data(raw_filename, 'header', hdr, 'chanindx', schan);
    
    stimonsample = find(diff(sdata)==24576)';
    stimoffsample = find(diff(sdata)==28672)';
    seqonsample = find(diff(sdata)==30720)';

 if subjectnr==12 % missing 2 triggers
        idx = find(diff(stimonsample)./hdr.Fs>.25 & diff(stimonsample)./hdr.Fs<.45);
        %[stimonsample(idx) stimonsample(idx+1)]./hdr.Fs;
        stimonsample = sort([stimonsample; stimonsample(idx)+mode(diff(stimonsample))]);
    end
% I think below is to make sure the photodiode # of events matches the
% expected number
    assert(numel(stimonsample)<=size(T,1),'found too many events!')
    assert(numel(stimonsample)>=size(T,1),'found not enough events!')


 %% changing the decoding time period 

    prestim = 0.1; %time period before stimulus is presented (in seconds)
    poststim = 0.8; %time period after stimulus is presented (in seconds)
    TRL = [stimonsample-ceil(prestim*hdr.Fs),ceil(stimonsample+poststim*hdr.Fs)]; %ciel rounds to the nearest integer
    TRL(:,3) = TRL(:,1)-stimonsample;

%% behavioral information 

    T = readtable(source_behavfn);
    T.subjectnr = T.eventnumber*0+subjectnr; %adds subject number to the farthest column 
    T2=table();
    T2.onset = stimonsample./hdr.Fs; %time of onset
    T2.duration(:) = floor(.100*hdr.Fs); %duration of the trial - floor rounds down
    T2.onsetsample = stimonsample; % sample number corresponding to the time of onset
    T = [T2 T] %combines T and T2 to T

    xd = diff(diff([T.onset T.time_stimon])')';
    assert(all(xd<.11 | xd>1),'event times do not seem to match')
    
    writetable(T,behavfn,'Filetype','text','delimiter','\t')
    
%% preprocessing 
    layout = ft_prepare_layout(struct('layout','biosemi64.lay')); % tells it what EEG system data was collected on 
    cfg = [];
    cfg.datafile = raw_filename; %specifies the name of the raw file
    cfg.trl = TRL; % from TRL table 
    %cfg.trl = TRL(1:10,:);
    %cfg.channel = 1:64;
    cfg.channel = hdr.label(1:64); %labels of the EEG channels
    cfg.reref = 'yes';
    cfg.refchannel = hdr.label(1:64);
    cfg.lpfilter = 'yes'; % want low-pass filter
    cfg.lpfreq = 100; % frequency for the low-pass filter (CHANGED FROM LPFILTER - we didn't actually filter our data in the paper)
    cfg.hpfilter = 'yes'; % want high-pass filter
    cfg.hpfreq = .1; % frequency for high-pass filter (CHANGED FROM HPFILTER - we didn't actually filter our data in the paper)
    cfg.hpfiltord = 4;
    cfg.dftfilter = 'yes'; % want to enable notch filtering to eliminate power line noise
    data_raw = ft_preprocessing(cfg); %the actual command to start preprocessing of cfg dataset with everything specified above (trials, filtering)
    data_raw.orig_label = data_raw.label; %labels data with original channels
    data_raw.label(1:64) = layout.label(1:64); %not sure why this is needed/how it is different from line above

%% resampling 
    cfg=[];
    cfg.resamplefs = 200; %frequency we are downsampling to 
    cfg.detrend = 'no'; % we are not detrending
    cfg.demean = 'yes'; % we are demeaning the data (baseline correction)
    cfg.baselinewindow = [-0.1 0]; %window that the data is demeaned to (in seconds - so it's demeaning to the 100 ms before stimulus) 
    data = ft_resampledata(cfg,data_raw); %command that actually resamples the data according the specifications above

%% cosmo
    cfg = [];
    cfg.keeptrials = 'yes'; % so individual trials are returned instead of averaged together
    ft = ft_timelockanalysis(cfg,data); % timelockes trials and creates covariance structure for each trial
    ds = cosmo_meeg_dataset(ft); %returns a dataset structure from timelocked pre-processed data
    ds.sa = table2struct(T,'toscalar',1); %sa = sample attributes (need to define this) 
    cosmo_check_dataset(ds); %checks dataset to make sure nothing is missing
    
%% save
    fprintf('saving...\n') %tells you it's saving
    save(cosmofn,'ds','-v7.3') %saving the dataset! 
    fprintf('done\n') %tells you when it's done with the code







