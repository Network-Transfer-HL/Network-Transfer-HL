% Create Timing Files for Delayed Recall Memory Analysis
%
% Author:   Alexandra Sobczak, M.Sc.
% Email:    alexandra.sobczak@uni-luebeck.de
% Date:     2023-03-14 (YYYY/MM/DD)
% Institute:University of Luebeck, IPSY1, Bunzeck Lab
% Project:  NetzTran

clear all
close all
clc

workpath_fMRI = ('/...'); % Pfad zu den fMRT Daten
workpath_behav = ('/...'); % Pfad zu den Daten von den Lernaufgaben
workpath_log = ('/...'); % Pfad zu den Logfiles

% We want to analyze the cue phase of the stimuli (figure pairs and word pairs) during delayed recall

% 2 tasks (FPA and NSWP) with 2 runs per task -> 4 runs in total
% 120 trials in total, 60 per task; 30 trials per run
% FPA and NSWP will be coded separately
% all events will be modeled: cue, target, feedback
% responses are coded as follows:
% 0: falsche Antwort
% 1: richtige Antwort
% 5: doppelt (also cue oder target aus Encode oder IR oder DR ist doppelt)
% 8: nicht gelernt (also cue oder target aus IR oder DR sind nicht in Encode aufgetaucht)
% 9: nicht beantwortet
% -> in conclusion, we have 30 conditions: 2 (task) x 3 (event (cue, target, feedback)) x 5 (response code)
% these conditions will later appear in our design matrix in SPM

task = {'FPA' 'NSWP'};
event = {'cue' 'target' 'feedback'};

TR = 1.84; % every 1.84s (1840ms) a volume is acquired

% ---------------------------------------------------------- %
% all subjects
VPnr = {'2','3',...};
code = {'VP02','VP03',...}; % folder
CBBM_enco_imm = {'13927',...};
CBBM_del = {'13928','13958',...};

Task_A = {'F',...};
Task_B = {'N',...};


%%
% for K = 40
 for K = 28:size(code,2)
    
    info.VPnr = VPnr{K};
    
    % ------------------------------------------------------------------- %
    % ---------------------------- load files---------------------------- %
    % ------------------------------------------------------------------- %
    % load all the stuff we need: nifti info, excel files, and log-files
    
    if strcmp(Task_A{K},'F')==1
        info.order = 1;
        info.taskA = 'FPA';
        info.taskB = 'NSWP';
    elseif strcmp(Task_A{K},'N')==1
        info.order = 2;
        info.taskA = 'NSWP';
        info.taskB = 'FPA';
    end
    
    % ---------------------------- NIFTI FILES -------------------------- %
    % get number of images in nifti files
    % FPA
    cd([workpath_fMRI,'FPA/',CBBM_del{K},'/run1']) % go to subject folder
    tmp_1 = dir(fullfile('a*DelayedRecall_1*.nii'));
    filenames.nifti.FPA_run1 = tmp_1.name;
    enco_FPA_run1 = niftiinfo(filenames.nifti.FPA_run1);
    n_scans_nifti_FPA_run1 = enco_FPA_run1.ImageSize(4);
    
    cd([workpath_fMRI,'FPA/',CBBM_del{K},'/run2']) % go to subject folder
    tmp_2 = dir(fullfile('a*DelayedRecall_2*.nii'));
    filenames.nifti.FPA_run2 = tmp_2.name;
    enco_FPA_run2 = niftiinfo(filenames.nifti.FPA_run2);
    n_scans_nifti_FPA_run2 = enco_FPA_run2.ImageSize(4);
    
    % NSWP
    cd([workpath_fMRI,'NSWP/',CBBM_del{K},'/run1']) % go to subject folder
    tmp_3 = dir(fullfile('a*DelayedRecall_1*.nii'));
    filenames.nifti.NSWP_run1 = tmp_3.name;
    enco_NSWP_run1 = niftiinfo(filenames.nifti.NSWP_run1);
    n_scans_nifti_NSWP_run1 = enco_NSWP_run1.ImageSize(4);
    
    cd([workpath_fMRI,'NSWP/',CBBM_del{K},'/run2']) % go to subject folder
    tmp_4 = dir(fullfile('a*DelayedRecall_2*.nii'));
    filenames.nifti.NSWP_run2 = tmp_4.name;
    enco_NSWP_run2 = niftiinfo(filenames.nifti.NSWP_run2);
    n_scans_nifti_NSWP_run2 = enco_NSWP_run2.ImageSize(4);

    n_scans_nifti_FPA = n_scans_nifti_FPA_run1+n_scans_nifti_FPA_run2;
    n_scans_nifti_NSWP = n_scans_nifti_NSWP_run1+n_scans_nifti_NSWP_run2;    

    
    % ------------------------ BEHAVIORAL DATA -------------------------- %
    cd([workpath_behav,code{K}]) % go to path where behavioral data is stored
    
    % get name of excel files with behavioral responses for DELAYED recall
    filename = dir(fullfile('*FPA_delayedRecall.xlsx'));
    filenames_behav.FPA = filename.name;
    
    filename = dir(fullfile('*NSWP_delayedRecall.xlsx'));
    filenames_behav.NSWP = filename.name;
    
    % load behavioral data from DELAYED recall
    % FPA --> still need excel files with 0,1,5,8,9 coding from Charlotte
    output.responses.FPA_run1 = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet1');
    output.responses.FPA_run2 = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet2');
    
    [~,output.cue.FPA_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet1','A2:A31');
    [~,output.cue.FPA_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet2','A2:A31');
    
    [~,output.target.FPA_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet1','B2:B31');
    [~,output.target.FPA_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.FPA],'Sheet2','B2:B31');
    
    % NSWP
    output.responses.NSWP_run1 = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet1');
    output.responses.NSWP_run2 = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet2');
    
    [~,output.cue.NSWP_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet1','A2:A31');
    [~,output.cue.NSWP_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet2','A2:A31');
    
    [~,output.target.NSWP_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet1','B2:B31');
    [~,output.target.NSWP_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Sheet2','B2:B31');
% % %    Bis VP7 (K=5) "Tabelle1" und "Tabelle2", Danach "Sheet1" und "Sheet2"
%     % NSWP
%     output.responses.NSWP_run1 = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle1');
%     output.responses.NSWP_run2 = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle2');
%     
%     [~,output.cue.NSWP_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle1','A2:A31');
%     [~,output.cue.NSWP_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle2','A2:A31');
%     
%     [~,output.target.NSWP_run1] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle1','B2:B31');
%     [~,output.target.NSWP_run2] = xlsread([workpath_behav,code{K},'/',filenames_behav.NSWP],'Tabelle2','B2:B31');
%     
    % ---------------------------- LOG FILES ---------------------------- %
    cd([workpath_log,code{K}]) % go to path where lof-files are stored
    % get logfile names
    filename = dir(fullfile('*Abfrage_FPA*.log'));
    filename_log_FPA = filename.name;
    
    filename = dir(fullfile('*Abfrage_NSWP*.log'));
    filename_log_NSWP = filename.name;
    
    % load logfiles
    data_FPA = table2cell(readtable([workpath_log,code{K},'/',filename_log_FPA],'FileType','text'));
    data_NSWP = table2cell(readtable([workpath_log,code{K},'/',filename_log_NSWP],'FileType','text'));
    
    % define name of stimulus events in logfile
    stimCue_FPA = 'bottomRight: autoDraw = True'; % cue
    stimTarget_FPA = 'Pfeile: autoDraw = True'; % target response
    stimTargetOff_FPA = 'Pfeile: autoDraw = False'; % offset of target
    stimFeedback_FPA = 'Figure2: autoDraw = True'; % feedback
    
    stimCue_NSWP = 'bottomRightTxt: autoDraw = True'; % presentation
    stimTarget_NSWP = 'Pfeile: autoDraw = True'; % response
    stimTargetOff_NSWP = 'Pfeile: autoDraw = False'; % offset of target
    stimFeedback_NSWP = 'TargetWord: autoDraw = True'; % feedback
    
    clearvars filename*
    
    % ------------------------------------------------------------------ %
    % ----------------- get all relevant timings ----------------------- %
    % ------------------------------------------------------------------ %
    % The scanner was stopped between task A and task B as well as between
    % the two runs of each task. While the scanner was stopped, psychopy
    % was running continuously. In order to align the timing of the scanner
    % with the timings recorded by psychopy in the log-files, we need to
    % separate the log-file for the two runs from one another and align
    % each run to zero.
    % The number of recorded scans in the log-files does not match the
    % number of scans in the nifti-files. The block duration is determined
    % by the number of scans (*TR) in the nifti files. That means that extra scans at the end of the
    % log-file will be ignored and missing scans will be appended.
    % We want to concatenate all runs and act as if the scanner was continuously running.
    % That is, the next run will artifically start in succession to the previous block.
    
    % 1) Calculate run durations using TR and number of scans
    % I have double checked that the calculated timings based on TR match the
    % timings in the log-file (ataking into account first scan onset of a run
    % and number of scans, etc). This way, we can be sure that the timings of
    % the events/stimuli is correct in relation to the calculated scan onsets.
    dur_FPA_run1 = TR * n_scans_nifti_FPA_run1;
    dur_FPA_run2 = TR * n_scans_nifti_FPA_run2;
    dur_NSWP_run1 = TR * n_scans_nifti_NSWP_run1;
    dur_NSWP_run2 = TR * n_scans_nifti_NSWP_run2;
    
    % 2) Separate logfiles into runs
    wait_4_scan_FPA = find(strcmp(data_FPA(:,3),'Wait_for_Scanner: autoDraw = True'));
    wait_4_scan_NSWP = find(strcmp(data_NSWP(:,3),'Wait_for_Scanner: autoDraw = True'));
    
    data_FPA_run1 = data_FPA(wait_4_scan_FPA(1):wait_4_scan_FPA(2)-1,:);
    data_FPA_run2 = data_FPA(wait_4_scan_FPA(2):end,:);
    
    data_NSWP_run1 = data_NSWP(wait_4_scan_NSWP(1):wait_4_scan_NSWP(2)-1,:);
    data_NSWP_run2 = data_NSWP(wait_4_scan_NSWP(2):end,:);
    
    % 3) Separate runs into trials
    % + doppelte 'target off' events rauswerfen
    % FPA %
    % run1
    ind_new_trial = find(contains(data_FPA_run1(:,3),'New trial'));
    for i = 1:length(ind_new_trial)-1
        data_FPA_run1_trials{i} = data_FPA_run1(ind_new_trial(i):ind_new_trial(i+1)-1,:);
    end
    data_FPA_run1_trials{length(ind_new_trial)} = data_FPA_run1(ind_new_trial(end):end,:);
    
    start_FPA_run1 = data_FPA_run1(1:ind_new_trial(1)-1,:);
    
    clearvars ind_new_trial i
    
    % check count of 'target off' events in trial & delete doubling events if there are any
    for i = 1:length(data_FPA_run1_trials)
        ind_doubling = find(strcmp(data_FPA_run1_trials{i}(:,3),stimTargetOff_FPA)); % check count of 'target off' events in trial
        if length(ind_doubling) > 1 % if more than one  'target off' event in one trial
            data_FPA_run1_trials{i}(ind_doubling(2),:) = []; % delete the second 'target off' event from the trial
        end
    end
    
    % run2
    ind_new_trial = find(contains(data_FPA_run2(:,3),'New trial'));
    for i = 1:length(ind_new_trial)-1
        data_FPA_run2_trials{i} = data_FPA_run2(ind_new_trial(i):ind_new_trial(i+1)-1,:);
    end
    data_FPA_run2_trials{length(ind_new_trial)} = data_FPA_run2(ind_new_trial(end):end,:);
    
    start_FPA_run2 = data_FPA_run2(1:ind_new_trial(1)-1,:);
    
    clearvars ind_new_trial i
    
    % check count of 'target off' events in trial & delete doubling events if there are any
    for i = 1:length(data_FPA_run2_trials)
        ind_doubling = find(strcmp(data_FPA_run2_trials{i}(:,3),stimTargetOff_FPA)); % check count of 'target off' events in trial
        if length(ind_doubling) > 1 % if more than one  'target off' event in one trial
            data_FPA_run2_trials{i}(ind_doubling(2),:) = []; % delete the second 'target off' event from the trial
        end
    end
    
    % NSWP %
    % run1
    ind_new_trial = find(contains(data_NSWP_run1(:,3),'New trial'));
    for i = 1:length(ind_new_trial)-1
        data_NSWP_run1_trials{i} = data_NSWP_run1(ind_new_trial(i):ind_new_trial(i+1)-1,:);
    end
    data_NSWP_run1_trials{length(ind_new_trial)} = data_NSWP_run1(ind_new_trial(end):end,:);
    
    start_NSWP_run1 = data_NSWP_run1(1:ind_new_trial(1)-1,:);
    
    clearvars ind_new_trial i
    
%     % check count of 'target off' events in trial & delete doubling events if there are any
%     for i = 1:length(data_NSWP_run1_trials)
%         ind_doubling = find(strcmp(data_NSWP_run1_trials{i}(:,3),stimTargetOff_NSWP)); % check count of 'target off' events in trial
%         if length(ind_doubling) > 1 % if more than one  'target off' event in one trial
%             data_NSWP_run1_trials{i}(ind_doubling(2),:) = []; % delete the second 'target off' event from the trial
%         end
%     end
    % check count of 'target off' events in trial & delete doubling events if there are any
    for i = 1:length(data_NSWP_run1_trials)
        ind_doubling = find(strcmp(data_NSWP_run1_trials{i}(:,3), stimTargetOff_NSWP)); % check count of 'target off' events in trial
        if length(ind_doubling) > 1 % if more than one 'target off' event in one trial
        data_NSWP_run1_trials{i}(ind_doubling(2:end), :) = []; % delete the second and subsequent 'target off' events from the trial
        end
    end
    % run2
    ind_new_trial = find(contains(data_NSWP_run2(:,3),'New trial'));
    for i = 1:length(ind_new_trial)-1
        data_NSWP_run2_trials{i} = data_NSWP_run2(ind_new_trial(i):ind_new_trial(i+1)-1,:);
    end
    data_NSWP_run2_trials{length(ind_new_trial)} = data_NSWP_run2(ind_new_trial(end):end,:);
    
    start_NSWP_run2 = data_NSWP_run2(1:ind_new_trial(1)-1,:);
    
    clearvars ind_new_trial i
    
    % check count of 'target off' events in trial & delete doubling events if there are any
    for i = 1:length(data_NSWP_run2_trials)
        ind_doubling = find(strcmp(data_NSWP_run2_trials{i}(:,3),stimTargetOff_NSWP)); % check count of 'target off' events in trial
        if length(ind_doubling) > 1 % if more than one  'target off' event in one trial
            data_NSWP_run2_trials{i}(ind_doubling(2),:) = []; % delete the second 'target off' event from the trial
        end
    end
    
    % after elimination of doubling 'target off' events, we put the data back together
    % keep the old log file data (just in case)
    data_FPA_run1_old = data_FPA_run1;
    data_FPA_run2_old = data_FPA_run2;
    data_NSWP_run1_old = data_NSWP_run1;
    data_NSWP_run2_old = data_NSWP_run2;
    clearvars data_FPA_run1 data_FPA_run2 data_NSWP_run1 data_NSWP_run2
    
    % put data back together without doubling events 
    % and use this variable for the next steps
    data_FPA_run1 = [start_FPA_run1;vertcat(data_FPA_run1_trials{:})];
    data_FPA_run2 = [start_FPA_run2;vertcat(data_FPA_run2_trials{:})];
    data_NSWP_run1 = [start_NSWP_run1;vertcat(data_NSWP_run1_trials{:})];
    data_NSWP_run2 = [start_NSWP_run2;vertcat(data_NSWP_run2_trials{:})];
    
    
    % 4) Get onsets of scans from logfile
    Tscan_FPA_run1 = [data_FPA_run1{strcmp(data_FPA_run1(:,3),'Keypress: s'),1}];
    Tscan_FPA_run2 = [data_FPA_run2{strcmp(data_FPA_run2(:,3),'Keypress: s'),1}];
    
    Tscan_NSWP_run1 = [data_NSWP_run1{strcmp(data_NSWP_run1(:,3),'Keypress: s'),1}];
    Tscan_NSWP_run2 = [data_NSWP_run2{strcmp(data_NSWP_run2(:,3),'Keypress: s'),1}];
    
    
    % 5) Get timing of first scan of the run
    t0_FPA_run1 = Tscan_FPA_run1(1);
    t0_FPA_run2 = Tscan_FPA_run2(1);
    
    t0_NSWP_run1 = Tscan_NSWP_run1(1);
    t0_NSWP_run2 = Tscan_NSWP_run2(1);
    
    % 6) Align run data time stamps to first scan (t0) of each run
    % add a fourth column to store the aligned timings
    data_FPA_run1(:,4)  = num2cell(cell2mat(data_FPA_run1(:,1))-t0_FPA_run1);
    data_FPA_run2(:,4)  = num2cell(cell2mat(data_FPA_run2(:,1))-t0_FPA_run2);
    
    data_NSWP_run1(:,4) = num2cell(cell2mat(data_NSWP_run1(:,1))-t0_NSWP_run1);
    data_NSWP_run2(:,4) = num2cell(cell2mat(data_NSWP_run2(:,1))-t0_NSWP_run2);
    
    
    % 7) Get onsets of events/stimuli (cue, target, feedback)
    my_onsets.cue.FPA_run1 = [data_FPA_run1{strcmp(data_FPA_run1(:,3),stimCue_FPA),4}]';
    my_onsets.cue.FPA_run2 = [data_FPA_run2{strcmp(data_FPA_run2(:,3),stimCue_FPA),4}]';
    
    my_onsets.target.FPA_run1 = [data_FPA_run1{strcmp(data_FPA_run1(:,3),stimTarget_FPA),4}]';
    my_onsets.target.FPA_run2 = [data_FPA_run2{strcmp(data_FPA_run2(:,3),stimTarget_FPA),4}]';
    
    my_onsets.feedback.FPA_run1 = [data_FPA_run1{strcmp(data_FPA_run1(:,3),stimFeedback_FPA),4}]';
    my_onsets.feedback.FPA_run2 = [data_FPA_run2{strcmp(data_FPA_run2(:,3),stimFeedback_FPA),4}]';
    
    my_onsets.cue.NSWP_run1 = [data_NSWP_run1{strcmp(data_NSWP_run1(:,3),stimCue_NSWP),4}]';
    my_onsets.cue.NSWP_run2 = [data_NSWP_run2{strcmp(data_NSWP_run2(:,3),stimCue_NSWP),4}]';
    
    my_onsets.target.NSWP_run1 = [data_NSWP_run1{strcmp(data_NSWP_run1(:,3),stimTarget_NSWP),4}]';
    my_onsets.target.NSWP_run2 = [data_NSWP_run2{strcmp(data_NSWP_run2(:,3),stimTarget_NSWP),4}]';
    
    my_onsets.feedback.NSWP_run1 = [data_NSWP_run1{strcmp(data_NSWP_run1(:,3),stimFeedback_NSWP),4}]';
    my_onsets.feedback.NSWP_run2 = [data_NSWP_run2{strcmp(data_NSWP_run2(:,3),stimFeedback_NSWP),4}]';
    
    % 8) Get event durations
    % get offsets of target for calculation of duration
    my_onsets.targetOff.FPA_run1 = [data_FPA_run1{strcmp(data_FPA_run1(:,3),stimTargetOff_FPA),4}]';
    my_onsets.targetOff.FPA_run2 = [data_FPA_run2{strcmp(data_FPA_run2(:,3),stimTargetOff_FPA),4}]';
    
    my_onsets.targetOff.NSWP_run1 = [data_NSWP_run1{strcmp(data_NSWP_run1(:,3),stimTargetOff_NSWP),4}]';
    my_onsets.targetOff.NSWP_run2 = [data_NSWP_run2{strcmp(data_NSWP_run2(:,3),stimTargetOff_NSWP),4}]';
    
    % calculate target durations (in seconds)
    targetDur_FPA_run1 = my_onsets.targetOff.FPA_run1-my_onsets.target.FPA_run1;
    targetDur_FPA_run2 = my_onsets.targetOff.FPA_run2-my_onsets.target.FPA_run2;
    
    targetDur_NSWP_run1 = my_onsets.targetOff.NSWP_run1-my_onsets.target.NSWP_run1;
    targetDur_NSWP_run2 = my_onsets.targetOff.NSWP_run2-my_onsets.target.NSWP_run2;
    
    % cue, target, and feedback durations for all runs together
    n_trials = length(my_onsets.cue.FPA_run1)+length(my_onsets.cue.FPA_run2)+length(my_onsets.cue.NSWP_run1)+length(my_onsets.cue.NSWP_run2);
    cue_dur_all = repmat(3.0,n_trials,1); % cues were presented for 3s
    feedback_dur_all = repmat(3.0,n_trials,1); % feedback was presented for 3s
    
    % target was presented until response (4s max)
    % concatenate target duration according to the order of tasks
    if info.order == 1 % task A = FPA; task B = NSWP
        target_dur_all = [targetDur_FPA_run1; targetDur_FPA_run2; targetDur_NSWP_run1; targetDur_NSWP_run2];
    elseif info.order == 2 % task A = NSWP; task B = FPA
        target_dur_all = [targetDur_NSWP_run1; targetDur_NSWP_run2; targetDur_FPA_run1; targetDur_FPA_run2];
    end
    
    
    
    % ------------------------------------------------------------------ %
    % --------------------- concatenate run onsets --------------------- %
    % ------------------------------------------------------------------ %
    % Make timings continuous by adding the durations of the previous blocks to the onset times
    
    
    % concatenate scan onsets
    % add durations of previous runs to stimulus onset times
    run_dur_all = dur_FPA_run1 + dur_FPA_run2 + dur_NSWP_run1 + dur_NSWP_run2;
    
    if info.order == 1 % task A = FPA; task B = NSWP
        
        my_onsets.taskA = 'FPA';
        my_onsets.taskB = 'NSWP';
        
        for i = 1:size(event,2)
            
            tmp_run2 = my_onsets.(event{i}).FPA_run2+dur_FPA_run1;
            tmp_run3 = my_onsets.(event{i}).NSWP_run1+dur_FPA_run1+dur_FPA_run2;
            tmp_run4 = my_onsets.(event{i}).NSWP_run2+dur_FPA_run1+dur_FPA_run2+dur_NSWP_run1;
            
            my_onsets.(event{i}).concatenate = [my_onsets.(event{i}).FPA_run1;tmp_run2;tmp_run3;tmp_run4];
            
        end
        
    elseif info.order == 2 % task A = NSWP; task B = FPA
        
        my_onsets.taskA = 'NSWP';
        my_onsets.taskB = 'FPA';
        
        for i = 1:size(event,2)
            
            tmp_run2 = my_onsets.(event{i}).NSWP_run2+dur_NSWP_run1;
            tmp_run3 = my_onsets.(event{i}).FPA_run1+dur_NSWP_run1+dur_NSWP_run2;
            tmp_run4 = my_onsets.(event{i}).FPA_run2+dur_NSWP_run1+dur_NSWP_run2+dur_FPA_run1;
            
            my_onsets.(event{i}).concatenate = [my_onsets.(event{i}).NSWP_run1;tmp_run2;tmp_run3;tmp_run4];
            
        end
    end
    
    
    % ------------------------------------------------------------------ %
    % --------------------- assign condition labels -------------------- %
    % ------------------------------------------------------------------ %
    % assign condition labels based on response codes in 'output.responses'
    % concatenate FPA and NSWP in the right order
    % this array has the same order as the onsets in 'my_onsets.concatenate'
    condition_FPA = cell(length(my_onsets.cue.FPA_run1)+length(my_onsets.cue.FPA_run2),1);
    condition_NSWP = cell(length(my_onsets.cue.NSWP_run1)+length(my_onsets.cue.NSWP_run2),1);
    
    cond_tmp_FPA = [output.responses.FPA_run1;output.responses.FPA_run2];
    cond_tmp_NSWP = [output.responses.NSWP_run1;output.responses.NSWP_run2];
    
    condition_FPA(cond_tmp_FPA ==1,1)={'FPA correct'};
    condition_FPA(cond_tmp_FPA ==0,1)={'FPA incorrect'};
    condition_FPA(cond_tmp_FPA ==5,1)={'FPA invalid code 5'};
    condition_FPA(cond_tmp_FPA ==8,1)={'FPA invalid code 8'};
    condition_FPA(cond_tmp_FPA ==9,1)={'FPA invalid code 9'};
    
    condition_NSWP(cond_tmp_NSWP ==1,1)={'NSWP correct'};
    condition_NSWP(cond_tmp_NSWP ==0,1)={'NSWP incorrect'};
    condition_NSWP(cond_tmp_NSWP ==5,1)={'NSWP invalid code 5'};
    condition_NSWP(cond_tmp_NSWP ==8,1)={'NSWP invalid code 8'};
    condition_NSWP(cond_tmp_NSWP ==9,1)={'NSWP invalid code 9'};
    
    if info.order == 1 % task A = FPA; task B = NSWP
        condition = [condition_FPA;condition_NSWP];
    elseif info.order == 2 % task A = NSWP; task B = FPA
        condition = [condition_NSWP;condition_FPA];
    end
    
  % ------------------------------------------------------------------ %
    % ---------------------- create multicon-file ---------------------- %
    % ------------------------------------------------------------------ %
    % if you have multiple conditions, then entering the details one
    % condition at a time is very inefficient. The multicon option allows
    % to load all the required information in one go.
    
    % 30 conditions: 2 (task) x 3 (event) x 5 (response code)
    % FPA and NSWP will be coded separately
    % all events will be modeled: cue, target, feedback
    % responses are coded as follows:
    % 0: falsche Antwort
    % 1: richtige Antwort
    % 5: doppelt (also cue oder target aus Encode oder IR oder DR ist doppelt)
    % 8: nicht gelernt (also cue oder target aus IR oder DR sind nicht in Encode aufgetaucht)
    % 9: nicht beantwortet
    

    % do not model invalid trials
    names        = [];
    onsets       = [];
    durations    = [];
    
    % FPA - cue
    names{1}     = 'FPA cue correct';
    names{2}     = 'FPA cue incorrect';

    % FPA - target
    names{3}     = 'FPA target correct';
    names{4}     = 'FPA target incorrect';

%     % FPA - feedback
%     names{5}    = 'FPA feedback correct';
%     names{6}    = 'FPA feedback incorrect';

    % NSWP -cue
    names{7}    = 'NSWP cue correct';
    names{8}    = 'NSWP cue incorrect';

    % NSWP - target
    names{9}    = 'NSWP target correct';
    names{10}    = 'NSWP target incorrect';

    % NSWP - feedback
%     names{11}    = 'NSWP feedback correct';
%     names{12}    = 'NSWP feedback incorrect';

    
    
    % FPA - cue
    onsets{1}    = my_onsets.cue.concatenate(strcmp(condition,'FPA correct'));
    onsets{2}    = my_onsets.cue.concatenate(strcmp(condition,'FPA incorrect'));

    % FPA - target
    onsets{3}    = my_onsets.target.concatenate(strcmp(condition,'FPA correct'));
    onsets{4}    = my_onsets.target.concatenate(strcmp(condition,'FPA incorrect'));

%     % FPA - feedback
%     onsets{5}    = my_onsets.feedback.concatenate(strcmp(condition,'FPA correct'));
%     onsets{6}    = my_onsets.feedback.concatenate(strcmp(condition,'FPA incorrect'));

    % NSWP - cue
    onsets{7}    = my_onsets.cue.concatenate(strcmp(condition,'NSWP correct'));
    onsets{8}    = my_onsets.cue.concatenate(strcmp(condition,'NSWP incorrect'));

    % NSWP - target
    onsets{9}    = my_onsets.target.concatenate(strcmp(condition,'NSWP correct'));
    onsets{10}    = my_onsets.target.concatenate(strcmp(condition,'NSWP incorrect'));

%     % NSWP - feedback
%     onsets{11}    = my_onsets.feedback.concatenate(strcmp(condition,'NSWP correct'));
%     onsets{12}    = my_onsets.feedback.concatenate(strcmp(condition,'NSWP incorrect'));

    
    
    % FPA - cue
    durations{1}  = cue_dur_all(strcmp(condition,'FPA correct')); % cues were presented for 3s
    durations{2}  = cue_dur_all(strcmp(condition,'FPA incorrect'));

    % FPA - target
    durations{3}  = target_dur_all(strcmp(condition,'FPA correct')); % targets were presented until response
    durations{4}  = target_dur_all(strcmp(condition,'FPA incorrect'));

%     % FPA - feedback
%     durations{5}  = feedback_dur_all(strcmp(condition,'FPA correct')); % feedback was presented for 3s
%     durations{6}  = feedback_dur_all(strcmp(condition,'FPA incorrect'));

    % NSWP - cue
    durations{7}  = cue_dur_all(strcmp(condition,'NSWP correct')); % cues were presented for 3s
    durations{8}  = cue_dur_all(strcmp(condition,'NSWP incorrect'));

    % NSWP - target
    durations{9}  = target_dur_all(strcmp(condition,'NSWP correct')); % targets were presented until response
    durations{10}  = target_dur_all(strcmp(condition,'NSWP incorrect'));

%     % NSWP - feedback
%     durations{11}  = feedback_dur_all(strcmp(condition,'NSWP correct')); % feedback was presented for 3s
%     durations{12}  = feedback_dur_all(strcmp(condition,'NSWP incorrect'));


    save([workpath_behav,code{K},'/',VPnr{K},'_concatenated_onsets_MEMORY_delayed_recall_FPA_NSWP_combined_multicon-file.mat'],'names','onsets','durations');
    
    
    
    
    ONSETS_CUE_FPA_NSWP_combined = table(my_onsets.cue.concatenate, cue_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
    ONSETS_TARGET_FPA_NSWP_combined = table(my_onsets.target.concatenate, target_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
%     ONSETS_FEEDBACK_FPA_NSWP_combined = table(my_onsets.feedback.concatenate, feedback_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
    
    save([workpath_behav,code{K},'/',VPnr{K},'_concatenated_onsets_MEMORY_delayed_recall_FPA_NSWP_combined.mat'],'ONSETS_CUE_FPA_NSWP_combined','ONSETS_TARGET_FPA_NSWP_combined');
%     save([workpath_behav,code{K},'/',VPnr{K},'_concatenated_onsets_MEMORY_delayed_recall_FPA_NSWP_combined.mat'],'ONSETS_CUE_FPA_NSWP_combined','ONSETS_TARGET_FPA_NSWP_combined','ONSETS_FEEDBACK_FPA_NSWP_combined');
    
    clearvars -except K workpath_fMRI workpath_behav workpath_log task event TR VPnr code CBBM_enco_imm CBBM_del Task_A Task_B Enco_version ImmRec_version DelRec_version 
     
end

