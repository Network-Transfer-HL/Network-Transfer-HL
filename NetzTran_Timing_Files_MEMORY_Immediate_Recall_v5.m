% Create Timing Files for Immediate Recall Memory Analysis
%
% Author:   Alexandra Sobczak, M.Sc.
% Email:    alexandra.sobczak@uni-luebeck.de
% Date:     2023-03-14 (YYYY/MM/DD)
% Institute:University of Luebeck, IPSY1, Bunzeck Lab
% Project:  NetzTran

clear all
close all
clc

% workpath_fMRI = ('/Users/AlexandraSobczak/Documents/Projects/NetzTran/1st Level Analysis/ImmRecall/'); % Pfad von Alex
workpath_fMRI = ('/Volumes/mehr platz/NetzTran/fMRI/1st_level_analysis/ImmRecall/'); % Pfad von Charlotte
% workpath_behav = '/Users/AlexandraSobczak/Documents/Projects/NetzTran/'; % Pfad von Alex
workpath_behav = ('/Volumes/mehr platz/NetzTran/Behaviour/Automatisierung_Onsets_DM/'); % Pfad von Charlotte
% workpath_log = '/Users/AlexandraSobczak/Documents/Projects/NetzTran/'; % Pfad von Alex
workpath_log = ('/Volumes/mehr platz/NetzTran/Behaviour/'); % Charlottes Pfad

% We want to analysze the cue phase of the stimuli (figure pairs and word pairs) during immediate recall

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
% %for example subject
% % dummy is VP46; files have been copied to a separate folder
% VPnr = {'dummy'}; % % hier müssen noch alle Vp-Nummern ergänzt werden 1 bis 53 (exkl. der Probanden, die von der Analyse ausgeschlossen werden sollen)
% code = {'dummy_F'}; % hier müssen noch alle Ordner-Namen ergänzt werden (exkl. der Probanden, die von der Analyse ausgeschlossen werden sollen)
% CBBM_enco_imm = {'14452'};
% CBBM_del = {'14453'};
% Task_A = {'F','F'};
% Task_B = {'N','N'};
% % these information are only stored here for the sake of documentation and is not used for anything (and therefore commented out)
% % Enco_version = {6, 6}; % only for information; this information is not used for anything
% % ImmRec_version = {6, 6}; % only for information; this information is not used for anything
% % DelRec_version = {6, 6}; % only for information; this information is not used for anything

% ---------------------------------------------------------- %
% % subjects 12-end
% VPnr = {'12','13','14','15','16','18','19','20','21','22','24','25','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
% code = {'VP12','VP13','VP14','VP15','VP16','VP18','VP19','VP20','VP21','VP22','VP24','VP25','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
% CBBM_enco_imm = {'14111','14131','14141','14221','14188','14225','14239','14259','14265','14278','14312','14287','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
% CBBM_del = {'14112','14132','14142','14222','14189','14226','14240','14260','14266','14279','14313','14288','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};
% 
% Task_A = {'N','F','N','F','N','N','F','F','N','N','N','F','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
% Task_B = {'F','N','F','N','F','F','N','N','F','F','F','N','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};

% ---------------------------------------------------------- %
% all subjects
VPnr = {'2','3','4','6','7','8','9','11','12','13','14','15','16','18','19','20','21','22','24','25','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
code = {'VP02','VP03','VP04','VP06','VP07','VP08','VP09','VP11','VP12','VP13','VP14','VP15','VP16','VP18','VP19','VP20','VP21','VP22','VP24','VP25','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
CBBM_enco_imm = {'13927','13957','13970','14035','14066','14074','14088','14105','14111','14131','14141','14221','14188','14225','14239','14259','14265','14278','14312','14287','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
CBBM_del = {'13928','13958','13971','14036','14067','14075','14089','14106','14112','14132','14142','14222','14189','14226','14240','14260','14266','14279','14313','14288','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};

Task_A = {'F','F','F','N','N','N','F','F','N','F','N','F','N','N','F','F','N','N','N','F','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
Task_B = {'N','N','N','F','F','F','N','N','F','N','F','N','F','F','N','N','F','F','F','N','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};

VPnr = {'2','3','4','6','7','8','9','11','12','13','14','15','18','19','20','21','22','24','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
code = {'VP02','VP03','VP04','VP06','VP07','VP08','VP09','VP11','VP12','VP13','VP14','VP15','VP18','VP19','VP20','VP21','VP22','VP24','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
CBBM_enco_imm = {'13927','13957','13970','14035','14066','14074','14088','14105','14111','14131','14141','14221','14225','14239','14259','14265','14278','14312','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
CBBM_del = {'13928','13958','13971','14036','14067','14075','14089','14106','14112','14132','14142','14222','14226','14240','14260','14266','14279','14313','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};

Task_A = {'F','F','F','N','N','N','F','F','N','F','N','F','N','F','F','N','N','N','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
Task_B = {'N','N','N','F','F','F','N','N','F','N','F','N','F','N','N','F','F','F','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};

% % ---------------------------------------------------------- %
% % complete and useable, 
% % all except
% % 16 (K=10)
% % 25 (K=17): CBBM 14287 missing
% % 
% VPnr = {'2','3','4','8','9','12','13','14','15','16','18','19','20','21','22','24','25','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
% code = {'VP02','VP03','VP04','VP08','VP09','VP12','VP13','VP14','VP15','VP16','VP18','VP19','VP20','VP21','VP22','VP24','VP25','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
% CBBM_enco_imm = {'13927','13957','13970','14074','14088','14111','14131','14141','14221','14188','14225','14239','14259','14265','14278','14312','14287','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
% CBBM_del = {'13928','13958','13971','14075','14089','14112','14132','14142','14222','14189','14226','14240','14260','14266','14279','14313','14288','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};
% 
% Task_A = {'F','F','F','N','F','N','F','N','F','N','N','F','F','N','N','N','F','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
% Task_B = {'N','N','N','F','N','F','N','F','N','F','F','N','N','F','F','F','N','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};
% 
% % Versions of the tests need to be corrected for some VP (see overview and supplementary infos Charlotte)
% % these information are only stored here for the sake of documentation and is not used for anything (and therefore commented out)
% % % Enco_version = {1,1,1,1,1,2,5,3,4,3,2,4,2,2,5,6,6,4,3,4,5,5,4,4,5,3,3,3,5,4,6,2,5,4,6,6,5,2,3,3,2,6};
% % % ImmRec_version = {1,1,1,1,1,2,5,3,4,3,2,4,2,2,5,6,6,4,3,4,5,5,4,4,5,3,3,3,5,4,6,2,5,4,6,6,5,2,3,3,2,6};
% % % DelRec_version = {1,1,1,1,1,2,5,3,4,3,2,4,2,2,5,6,6,4,3,4,5,5,4,4,5,3,3,3,5,4,6,2,5,4,6,6,5,2,3,3,2,6};
% % ---------------------------------------------------------- %

%%
% for K = 17
 for K = 6:size(code,2)
    
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
    cd([workpath_fMRI,'FPA/',CBBM_enco_imm{K},'/run1']) % go to subject folder
    tmp_1 = dir(fullfile('a*ImmRecall_1*.nii'));
    filenames.nifti.FPA_run1 = tmp_1.name;
    enco_FPA_run1 = niftiinfo(filenames.nifti.FPA_run1);
    n_scans_nifti_FPA_run1 = enco_FPA_run1.ImageSize(4);
    
    cd([workpath_fMRI,'FPA/',CBBM_enco_imm{K},'/run2']) % go to subject folder
    tmp_2 = dir(fullfile('a*ImmRecall_2*.nii'));
    filenames.nifti.FPA_run2 = tmp_2.name;
    enco_FPA_run2 = niftiinfo(filenames.nifti.FPA_run2);
    n_scans_nifti_FPA_run2 = enco_FPA_run2.ImageSize(4);
    
    % NSWP
    cd([workpath_fMRI,'NSWP/',CBBM_enco_imm{K},'/run1']) % go to subject folder
    tmp_3 = dir(fullfile('a*ImmRecall_1*.nii'));
    filenames.nifti.NSWP_run1 = tmp_3.name;
    enco_NSWP_run1 = niftiinfo(filenames.nifti.NSWP_run1);
    n_scans_nifti_NSWP_run1 = enco_NSWP_run1.ImageSize(4);
    
    cd([workpath_fMRI,'NSWP/',CBBM_enco_imm{K},'/run2']) % go to subject folder
    tmp_4 = dir(fullfile('a*ImmRecall_2*.nii'));
    filenames.nifti.NSWP_run2 = tmp_4.name;
    enco_NSWP_run2 = niftiinfo(filenames.nifti.NSWP_run2);
    n_scans_nifti_NSWP_run2 = enco_NSWP_run2.ImageSize(4);

    n_scans_nifti_FPA = n_scans_nifti_FPA_run1+n_scans_nifti_FPA_run2;
    n_scans_nifti_NSWP = n_scans_nifti_NSWP_run1+n_scans_nifti_NSWP_run2;    

    
    % ------------------------ BEHAVIORAL DATA -------------------------- %
    cd([workpath_behav,code{K}]) % go to path where behavioral data is stored
    
    % get name of excel files with behavioral responses for IMMEDIATE recall
    filename = dir(fullfile('*FPA_immediateRecall.xlsx'));
    filenames_behav.FPA = filename.name;
    
    filename = dir(fullfile('*NSWP_immediateRecall.xlsx'));
    filenames_behav.NSWP = filename.name;
    
    % load behavioral data from IMMEDIATE recall
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
    
    % ---------------------------- LOG FILES ---------------------------- %
    cd([workpath_log,code{K}]) % go to path where lof-files are stored
    % get logfile names
    filename = dir(fullfile('*ImmediateRecall_FPA*.log'));
    filename_log_FPA = filename.name;
    
    filename = dir(fullfile('*ImmediateRecall_NSWP*.log'));
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
    
    % check count of 'target off' events in trial & delete doubling events if there are any
    for i = 1:length(data_NSWP_run1_trials)
        ind_doubling = find(strcmp(data_NSWP_run1_trials{i}(:,3),stimTargetOff_NSWP)); % check count of 'target off' events in trial
        if length(ind_doubling) > 1 % if more than one  'target off' event in one trial
            data_NSWP_run1_trials{i}(ind_doubling(2),:) = []; % delete the second 'target off' event from the trial
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

    % FPA - feedback
    names{5}    = 'FPA feedback correct';
    names{6}    = 'FPA feedback incorrect';

    % NSWP -cue
    names{7}    = 'NSWP cue correct';
    names{8}    = 'NSWP cue incorrect';

    % NSWP - target
    names{9}    = 'NSWP target correct';
    names{10}    = 'NSWP target incorrect';

    % NSWP - feedback
    names{11}    = 'NSWP feedback correct';
    names{12}    = 'NSWP feedback incorrect';

    
    
    % FPA - cue
    onsets{1}    = my_onsets.cue.concatenate(strcmp(condition,'FPA correct'));
    onsets{2}    = my_onsets.cue.concatenate(strcmp(condition,'FPA incorrect'));

    % FPA - target
    onsets{3}    = my_onsets.target.concatenate(strcmp(condition,'FPA correct'));
    onsets{4}    = my_onsets.target.concatenate(strcmp(condition,'FPA incorrect'));

    % FPA - feedback
    onsets{5}    = my_onsets.feedback.concatenate(strcmp(condition,'FPA correct'));
    onsets{6}    = my_onsets.feedback.concatenate(strcmp(condition,'FPA incorrect'));

    % NSWP - cue
    onsets{7}    = my_onsets.cue.concatenate(strcmp(condition,'NSWP correct'));
    onsets{8}    = my_onsets.cue.concatenate(strcmp(condition,'NSWP incorrect'));

    % NSWP - target
    onsets{9}    = my_onsets.target.concatenate(strcmp(condition,'NSWP correct'));
    onsets{10}    = my_onsets.target.concatenate(strcmp(condition,'NSWP incorrect'));

    % NSWP - feedback
    onsets{11}    = my_onsets.feedback.concatenate(strcmp(condition,'NSWP correct'));
    onsets{12}    = my_onsets.feedback.concatenate(strcmp(condition,'NSWP incorrect'));

    
    
    % FPA - cue
    durations{1}  = cue_dur_all(strcmp(condition,'FPA correct')); % cues were presented for 3s
    durations{2}  = cue_dur_all(strcmp(condition,'FPA incorrect'));

    % FPA - target
    durations{3}  = target_dur_all(strcmp(condition,'FPA correct')); % targets were presented until response
    durations{4}  = target_dur_all(strcmp(condition,'FPA incorrect'));

    % FPA - feedback
    durations{5}  = feedback_dur_all(strcmp(condition,'FPA correct')); % feedback was presented for 3s
    durations{6}  = feedback_dur_all(strcmp(condition,'FPA incorrect'));

    % NSWP - cue
    durations{7}  = cue_dur_all(strcmp(condition,'NSWP correct')); % cues were presented for 3s
    durations{8}  = cue_dur_all(strcmp(condition,'NSWP incorrect'));

    % NSWP - target
    durations{9}  = target_dur_all(strcmp(condition,'NSWP correct')); % targets were presented until response
    durations{10}  = target_dur_all(strcmp(condition,'NSWP incorrect'));

    % NSWP - feedback
    durations{11}  = feedback_dur_all(strcmp(condition,'NSWP correct')); % feedback was presented for 3s
    durations{12}  = feedback_dur_all(strcmp(condition,'NSWP incorrect'));


    save([workpath_behav,code{K},'/',VPnr{K},'_concatenated_onsets_MEMORY_immediate_recall_FPA_NSWP_combined_multicon-file.mat'],'names','onsets','durations');
    
    
    
    
    ONSETS_CUE_FPA_NSWP_combined = table(my_onsets.cue.concatenate, cue_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
    ONSETS_TARGET_FPA_NSWP_combined = table(my_onsets.target.concatenate, target_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
    ONSETS_FEEDBACK_FPA_NSWP_combined = table(my_onsets.feedback.concatenate, feedback_dur_all, condition, 'VariableNames', { 'Onsets', 'Duration', 'Condition'} );
    
    save([workpath_behav,code{K},'/',VPnr{K},'_concatenated_onsets_MEMORY_immediate_recall_FPA_NSWP_combined.mat'],'ONSETS_CUE_FPA_NSWP_combined','ONSETS_TARGET_FPA_NSWP_combined','ONSETS_FEEDBACK_FPA_NSWP_combined');
    
    clearvars -except K workpath_fMRI workpath_behav workpath_log task event TR VPnr code CBBM_enco_imm CBBM_del Task_A Task_B Enco_version ImmRec_version DelRec_version 
    
end

