%%% Credits:
% % Author:   Alexandra Sobczak, M.Sc.
% Email:    alexandra.sobczak@uni-luebeck.de
% Date:     2023-03-14 (YYYY/MM/DD)
% Institute:University of Luebeck, IPSY1, Bunzeck Lab
% Project:  NetzTran
% % Co-Author/ Edited by: Charlotte Jeschina
% Email:   charlotte.jeschina@student.uni-luebeck.de
% Date:     2024

function statsBatch
% statistical analysis for fMRI data using SPM12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script 
% a) specifies the model,
% b) estimates the model
% c) contrasts
%
% See wiki on how to perform concatenation
% https://en.wikibooks.org/wiki/SPM/Concatenation
% 
% _________________________________________________________________________

%% Define subject parameters and directories
% ----------- Pfade definieren ---------- %
addpath /Applications/MATLAB_R2021b.app/toolbox/spm12

spm12dir = '/Applications/MATLAB_R2021b.app/toolbox/spm12';

fs          = filesep; %file sep
n_sess      = 1; % no of sessions (runs) this must be 1 for concatinated runs
workpath_fMRI = ('/Volumes/mehr platz/NetzTran/fMRI/1st_level_analysis/DelRecall/'); % Pfad von Charlotte
workpath_behav = ('/Volumes/mehr platz/NetzTran/Behaviour/computed_behaviour_data'); % Pfad von Charlotte
dir_base    = ('/Volumes/mehr platz/NetzTran/fMRI'); % Pfad von Charlotte
dir_results = '1st_Level/Memory_Delayed_Recall'; % Pfad zu den Ergebnissen der neuen 1st Level Analyse 
sess_prfx   = 'task';
TR = 1.84; % every 1.84s (1840ms) a volume is acquired

% % % % % ------------------------------------- %

% Weitere wichtige Pfade und Ordner erstellen
cd(dir_base)
mkdir 1st_Level
cd([dir_base fs '1st_Level'])
mkdir Memory_Delayed_Recall


% z.B. für VP 53
VPnr = {...,'53'};
code = {...,'VP53'};
CBBM_enco_imm = {...,'14612'};
CBBM_del = {...,'14613'};
Task_A = {...,'N'};
Task_B = {...,'F'};


spm fmri


%% Define what processing we want
specify  = 1;
concat   = 1;
estimate = 1;
contrast = 1;

%spm fmri
%% specify


if specify
    disp(['specify model ']);
    
    for s0 = 1 : length(code)
        disp(['Specify model for Subject: ', code{s0}]);
        
        % navigate to the directory for 1st Level results
        cd([dir_base fs dir_results]);
        
        % create a folder for the subject
        mkdir(code{s0})
        
        cd([dir_base fs dir_results fs code{s0} fs]);
        
        res_dir = pwd;
        b = cellstr([pwd]);
        jobs{1}.stats{1}.fmri_spec.dir = b;
        %set timing
        jobs{1}.stats{1}.fmri_spec.timing.units   = 'secs';
        jobs{1}.stats{1}.fmri_spec.timing.RT      = TR;
        jobs{1}.stats{1}.fmri_spec.timing.fmri_t  = 52; % vorher 28 (image size)
        jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 1;
        
        %set scans, cond, multi, regress, multi_reg, hrf
        %loop to define sessions in job
        for sess = 1:n_sess
            
            %%%%%% FPA %%%%%%
            %define epi's to be used in the session
            scanDir1 = [workpath_fMRI fs 'FPA/' CBBM_del{s0} fs 'run1' fs];
            %select scans and assign to job
            cd (scanDir1);
            f1   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelayedRecall_1.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_FPA_run1 = length(f1);
            
            scanDir2 = [workpath_fMRI fs 'FPA/' CBBM_del{s0} fs 'run2' fs];
            %select scans and assign to job
            cd (scanDir2);
            f2   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelayedRecall_2.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_FPA_run2 = length(f2);

            %%%%%% NSWP %%%%%%
            %define epi's to be used in the session
            scanDir3 = [workpath_fMRI fs 'NSWP/' CBBM_del{s0} fs 'run1' fs];
            %select scans and assign to job
            cd (scanDir3);
            f3   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelayedRecall_1.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_NSWP_run1 = length(f3);
            
            %define epi's to be used in the session
            scanDir4 = [workpath_fMRI fs 'NSWP/' CBBM_del{s0} fs 'run2' fs];
            %select scans and assign to job
            cd (scanDir4);
            f4   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelayedRecall_2.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_NSWP_run2 = length(f4);
            
            
            %check which task (FPA or NSWP) was performed first. 
            % The order of the files needs to match the order of tasks in the multicon-file!
            if strcmp(Task_A{s0},'F')==1 % if FPA is fist task, append files such that FPA files are first
                
                files  = [f1; f2; f3; f4];
            
            elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
                
                files  = [f3; f4; f1; f2];
                
            end
                
            jobs{1}.stats{1}.fmri_spec.sess(sess).scans = files;
            
            
            f = []; files = [];
            % load multicon_file
%             multicon_files = spm_select('List', [dir_base fs 'raw_data/MRI/young' fs code{s0} fs 'multicon_all-file.mat']);
%             multicon_file = cellstr([dir_base '/raw_data/MRI/young' fs code{s0} fs multicon_files(sess,:)]);

            multicon_files = spm_select('List', [workpath_behav fs code{s0} fs VPnr{s0} '_concatenated_onsets_MEMORY_delayed_recall_FPA_NSWP_combined_multicon-file.mat']);
            multicon_file = cellstr([workpath_behav fs code{s0} fs multicon_files(sess,:)]);
            
            
            jobs{1}.stats{1}.fmri_spec.sess(sess).multi = multicon_file;
            
            
            %load multi_reg (has to be created through concatination)
            %FPA 
            %-run1
            cd (scanDir1);
            tmp = spm_select('List', fullfile(['rp*DelayedRecall_1','*.txt']));
            a=load(tmp); clear tmp;
            %-run2
            cd (scanDir2);
            tmp = spm_select('List', fullfile(['rp*DelayedRecall_2','*.txt']));
            b=load(tmp); clear tmp;
            % NSWP
            %-run3
            cd (scanDir3);
            tmp = spm_select('List', fullfile(['rp*DelayedRecall_1','*.txt']));
            c=load(tmp); clear tmp;
            %-run4
            cd (scanDir4);
            tmp = spm_select('List', fullfile(['rp*DelayedRecall_2','*.txt']));
            d=load(tmp); clear tmp;

            %
            %check which task (FPA or NSWP) was performed first. 
            % The order of the files needs to match the order of tasks in the multicon-file!
            if strcmp(Task_A{s0},'F')==1 % if FPA is fist task, append files such that FPA files are first
                
                R  = [a; b; c; d];
            
            elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
                
                R  = [c; d; a; b];
                
            end
            
            
            
            cd([workpath_behav fs code{s0} fs]);
            
            
            save R;
            multi_reg = spm_select('List', pwd, 'R.mat');
            multi_reg = cellstr([pwd fs multi_reg]);
            jobs{1}.stats{1}.fmri_spec.sess(sess).multi_reg = multi_reg;
            jobs{1}.stats{1}.fmri_spec.sess(sess).hpf = 128;
            
        end
        
        jobs{1}.stats{1}.fmri_spec.fact = [];
        jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs = [1 1]; % uses hrf and  tmp derv and  time disp otherwise: if no [0 0]
        jobs{1}.stats{1}.fmri_spec.volt = 1;
        jobs{1}.stats{1}.fmri_spec.global = 'None';
        cvi = 'FAST'; % AR(1) not recommended for short TRs (<1.5s). Instead, use FAST.
        cellstr([cvi]);
        jobs{1}.stats{1}.fmri_spec.cvi = cvi;
        cd (res_dir);
        % save and run job
        save 1_specify_model.mat jobs
        disp(['SPECIFY MODEL for subject ' code{s0}]);
        spm_jobman('run','1_specify_model.mat');
        disp('Job completed');
        clear jobs;
        
    end
    
end



%% Concatinate runs: https://en.wikibooks.org/wiki/SPM/Concatenation
if concat
    
    % The order of the files needs to match the order of tasks in the multicon-file!
    if strcmp(Task_A{s0},'F')==1 % if FPA is fist task, append files such that FPA files are first
        
        scans = [n_scans_nifti_FPA_run1 n_scans_nifti_FPA_run2 n_scans_nifti_NSWP_run1 n_scans_nifti_NSWP_run2];
        
    elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
        
        scans = [n_scans_nifti_NSWP_run1 n_scans_nifti_NSWP_run2 n_scans_nifti_FPA_run1 n_scans_nifti_FPA_run2];
        
    end
    
    
    for s0 =  1:length(code)
        cd([dir_base fs dir_results fs code{s0}]);
        res_dir = pwd;
        f = spm_select('List', res_dir, 'SPM.mat');
        load(f);
        spm_fmri_concatenate('SPM.mat', scans)
        cd (res_dir);
    end
    
end


%% Estimate
if estimate
    disp(['estimate model ']);
    
    for s0 = 1 : length(code)
        disp(['Estimate model for Subject: ', code{s0}]);
        
        cd([dir_base fs dir_results fs code{s0}]);
        res_dir = pwd;
        %which directory
        f   = spm_select('List', res_dir, 'SPM.mat');
        spmmat = cellstr([res_dir fs f]);
        jobs{1}.stats{1}.fmri_est.spmmat = spmmat;
        
        %which method
        method = 'Classical: 1';
        method = cellstr([method]);
        jobs{1}.stats{1}.fmri_est.method = [1];
        cd (res_dir);
        
        % save and run job
        save 2_estimate_model.mat jobs;
        disp(['ESTIMATE MODEL for subject ' code{s0}]);
        spm_jobman('run','2_estimate_model.mat');
        disp('Job completed');
        clear jobs;
    end
end


%% Contrast
if contrast
    disp(['set contrasts ']);
    
    nregressors = 30 %3x FPA cue correct, 3x FPA cue incorrect, 3x FPA target correct, 3x FPA target incorrect, 3x NSWP cue correct, 3x NSWP cue incorrect, 3x NSWP target correct, 3x NSWP target incorrect, 6x Realignment
 
    
    for s0 = 1 : length(code)
        disp(['set contrasts for Subject: ', code{s0}]);
        
        cd([dir_base fs dir_results fs code{s0}]);
        res_dir = pwd;
        %which directory
        f   = spm_select('List', res_dir, 'SPM.mat');
        spmmat = cellstr([res_dir fs f]);
        jobs{1}.stats{1}.con.spmmat = spmmat;
        
        %define contrasts
        % The order of the files needs to match the order of tasks in the multicon-file!
        if strcmp(Task_A{s0},'F')==1 % if FPA is fist task, append files such that FPA files are first
            
            jobs{1}.stats{1}.con.consess{1}.tcon.name = 'FPA cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{1}.tcon.convec = repmat([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{1}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{2}.tcon.name = 'FPA cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{2}.tcon.convec = repmat([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{2}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{3}.tcon.name = 'FPA target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{3}.tcon.convec = repmat([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{3}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{4}.tcon.name = 'FPA target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{4}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{4}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{5}.tcon.name = 'NSWP cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{5}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{5}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{6}.tcon.name = 'NSWP cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{6}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{6}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{7}.tcon.name = 'NSWP target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{7}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{7}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{8}.tcon.name = 'NSWP target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{8}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{8}.tcon.sessrep = 'none';
            

            
            
        elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
            
            jobs{1}.stats{1}.con.consess{1}.tcon.name = 'FPA cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{1}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{1}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{2}.tcon.name = 'FPA cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{2}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{2}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{3}.tcon.name = 'FPA target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{3}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{3}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{4}.tcon.name = 'FPA target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{4}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{4}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{5}.tcon.name = 'NSWP cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{5}.tcon.convec = repmat([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{5}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{6}.tcon.name = 'NSWP cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{6}.tcon.convec = repmat([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{6}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{7}.tcon.name = 'NSWP target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{7}.tcon.convec = repmat([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{7}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{8}.tcon.name = 'NSWP target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{8}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{8}.tcon.sessrep = 'none';
            
            
        end
        
        jobs{1}.stats{1}.con.delete = 1; % delete existing contrast = 1
        cd (res_dir);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save and run job
        save 3_contrasts.mat jobs
        disp(['SETTING CONTRASTS for subject ' code{s0}]);
        spm_jobman('run','3_contrasts.mat');
        disp('Job completed');
        clear jobs
    end
end
