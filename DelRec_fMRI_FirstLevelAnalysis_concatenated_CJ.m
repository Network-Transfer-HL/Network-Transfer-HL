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
% addpath /Users/AlexandraSobczak/Documents/MATLAB/toolboxes/spm12 % Define SPM in MATLAB path
addpath /Applications/MATLAB_R2021b.app/toolbox/spm12

% spm12dir = '/Users/AlexandraSobczak/Documents/MATLAB/toolboxes/spm12';
spm12dir = '/Applications/MATLAB_R2021b.app/toolbox/spm12';

fs          = filesep; %file sep
n_sess      = 1; % no of sessions (runs) this must be 1 for concatinated runs
% workpath_fMRI = ('/Users/AlexandraSobczak/Documents/Projects/NetzTran/1st Level Analysis/Delecall'); % Pfad von Alex
workpath_fMRI = ('/Volumes/mehr platz/NetzTran/fMRI/1st_level_analysis/DelRecall/'); % Pfad von Charlotte
% workpath_behav = '/Users/AlexandraSobczak/Documents/Projects/NetzTran'; % Pfad von Alex
workpath_behav = ('/Volumes/mehr platz/NetzTran/Behaviour/computed_behaviour_data'); % Pfad von Charlotte
% dir_base    = '/Users/AlexandraSobczak/Documents/Projects/NetzTran'; % Pfad von Alex
dir_base    = ('/Volumes/mehr platz/NetzTran/fMRI'); % Pfad von Charlotte
dir_results = '1st_Level/Memory_Delayed_Recall'; % Pfad zu den Ergebnissen der neuen 1st Level Analyse (bei uns beiden gleich)
sess_prfx   = 'task';
TR = 1.84; % every 1.84s (1840ms) a volume is acquired

% % % % % ------------------------------------- %


% Weitere wichtige Pfade und Ordner erstellen
cd(dir_base)
mkdir 1st_Level
cd([dir_base fs '1st_Level'])
mkdir Memory_Immediate_Recall
mkdir Memory_Delayed_Recall
mkdir DM_Immediate_Recall
mkdir DM_Delayed_Recall



% 
%for example subject
% dummy is VP46; files have been copied to a separate folder
VPnr = {'26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
code = {'VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
CBBM_enco_imm = {'14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
CBBM_del = {'14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};
Task_A = {'F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
Task_B = {'N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};

% % Script hat gestoppt bei 16,25
% VPnr = {'2','3','4','6','7','8','9','11','12','13','14','15','16','18','19','20','21','22','24','25',
% code = {'VP02','VP03','VP04','VP06','VP07','VP08','VP09','VP11','VP12','VP13','VP14','VP15','VP16','VP18','VP19','VP20','VP21','VP22','VP24','VP25',
% CBBM_enco_imm = {'13927','13957','13970','14035','14066','14074','14088','14105','14111','14131','14141','14221','14188','14225','14239','14259','14265','14278','14312','14287',
% CBBM_del = {'13928','13958','13971','14036','14067','14075','14089','14106','14112','14132','14142','14222','14189','14226','14240','14260','14266','14279','14313','14288',
%  Task_A = {'F','F','F','N','N','N','F','F','N','F','N','F','N','N','F','F','N','N','N','F',
%  Task_B = {'N','N','N','F','F','F','N','N','F','N','F','N','F','F','N','N','F','F','F','N',


% VPnr = {'2','3','4','6','7','8','9','11','12','13','14','15','18','19','20','21','22','24','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
% code = {'VP02','VP03','VP04','VP06','VP07','VP08','VP09','VP11','VP12','VP13','VP14','VP15','VP18','VP19','VP20','VP21','VP22','VP24','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'};
% CBBM_enco_imm = {'13927','13957','13970','14035','14066','14074','14088','14105','14111','14131','14141','14221','14225','14239','14259','14265','14278','14312','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
% CBBM_del = {'13928','13958','13971','14036','14067','14075','14089','14106','14112','14132','14142','14222','14226','14240','14260','14266','14279','14313','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};
% Task_A = {'F','F','F','N','N','N','F','F','N','F','F','N','N','F','F','N','N','N','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
% Task_B = {'N','N','N','F','F','F','N','N','F','N','N','F','F','N','N','F','F','F','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};


% VPnr = {'12','13','14','15','16','18','19','20','21','22','24','25','26','27','28','30','31','32','33','34','38','40','41','42','43','44','45','46','47','48','49','50','52','53'};
% code = {'VP12','VP13','VP14','VP15','VP16','VP18','VP19','VP20','VP21','VP22','VP24','VP25','VP26','VP27','VP28','VP30','VP31','VP32','VP33','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP48','VP49','VP50','VP52','VP53'}; % folder
% CBBM_enco_imm = {'14111','14131','14141','14221','14188','14225','14239','14259','14265','14278','14312','14287','14338','14343','14345','14365','14363','14373','14388','14381','14409','14428','14421','14425','14466','14437','14482','14452','14479','14505','14519','14509','14542','14612'};
% CBBM_del = {'14112','14132','14142','14222','14189','14226','14240','14260','14266','14279','14313','14288','14339','14344','14346','14366','14364','14374','14389','14382','14410','14429','14422','14426','14467','14438','14483','14453','14480','14506','14520','14510','14543','14613'};
%  Task_A = {'N','F','N','F','N','N','F','F','N','N','N','F','F','F','N','N','N','N','N','N','N','F','N','F','N','N','N','F','F','F','F','N','F','N'};
%  Task_B = {'F','N','F','N','F','F','N','N','F','F','F','N','N','N','F','F','F','F','F','F','F','N','F','N','F','F','F','N','N','N','N','F','N','F'};


% VPnr = {'dummy'}; % % hier müssen noch alle Vp-Nummern ergänzt werden 1 bis 53 (exkl. der Probanden, die von der Analyse ausgeschlossen werden sollen)
% code = {'dummy_F'}; % hier müssen noch alle Ordner-Namen ergänzt werden (exkl. der Probanden, die von der Analyse ausgeschlossen werden sollen)
% CBBM_enco_imm = {'14452'};
% CBBM_del = {'14453'};
% Task_A = {'F','F'};
% Task_B = {'N','N'};
% these information are only stored here for the sake of documentation and is not used for anything (and therefore commented out)
% Enco_version = {6, 6}; % only for information; this information is not used for anything
% ImmRec_version = {6, 6}; % only for information; this information is not used for anything
% DelRec_version = {6, 6}; % only for information; this information is not used for anything


spm fmri

%% Define what processing we want
specify  = 1;
concat   = 1;
estimate = 1;
contrast = 1;

%spm fmri
%% specify
%
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
            f1   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelRecall_1.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_FPA_run1 = length(f1);
            
            scanDir2 = [workpath_fMRI fs 'FPA/' CBBM_del{s0} fs 'run2' fs];
            %select scans and assign to job
            cd (scanDir2);
            f2   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelRecall_2.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_FPA_run2 = length(f2);

            %%%%%% NSWP %%%%%%
            %define epi's to be used in the session
            scanDir3 = [workpath_fMRI fs 'NSWP/' CBBM_del{s0} fs 'run1' fs];
            %select scans and assign to job
            cd (scanDir3);
            f3   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelRecall_1.*.nii', inf)); %goes through all subfolders
            n_scans_nifti_NSWP_run1 = length(f3);
            
            %define epi's to be used in the session
            scanDir4 = [workpath_fMRI fs 'NSWP/' CBBM_del{s0} fs 'run2' fs];
            %select scans and assign to job
            cd (scanDir4);
            f4   = cellstr(spm_select('ExtFPListRec', pwd, '^swa.*.DelRecall_2.*.nii', inf)); %goes through all subfolders
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
            tmp = spm_select('List', fullfile(['rp*DelRecall_1','*.txt']));
            a=load(tmp); clear tmp;
            %-run2
            cd (scanDir2);
            tmp = spm_select('List', fullfile(['rp*DelRecall_2','*.txt']));
            b=load(tmp); clear tmp;
            % NSWP
            %-run3
            cd (scanDir3);
            tmp = spm_select('List', fullfile(['rp*DelRecall_1','*.txt']));
            c=load(tmp); clear tmp;
            %-run4
            cd (scanDir4);
            tmp = spm_select('List', fullfile(['rp*DelRecall_2','*.txt']));
            d=load(tmp); clear tmp;

            %
            %check which task (FPA or NSWP) was performed first. 
            % The order of the files needs to match the order of tasks in the multicon-file!
            if strcmp(Task_A{s0},'F')==1 % if FPA is fist task, append files such that FPA files are first
                
                R  = [a; b; c; d];
            
            elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
                
                R  = [c; d; a; b];
                
            end
            
            
%             cd([dir_base '/raw_data/MRI/young' fs code{s0} fs]);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    nregressors = 42 %3x FPA cue correct, 3x FPA cue incorrect, 3x FPA target correct, 3x FPA target incorrect, 3x FPA feedback correct, 3x FPA feedback incorrect, 3x NSWP cue correct, 3x NSWP cue incorrect, 3x NSWP target correct, 3x NSWP target incorrect, 3x NSWP feedback correct, 3x NSWP feedback incorrect, 6x Realignment
 
    
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
            jobs{1}.stats{1}.con.consess{1}.tcon.convec = repmat([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{1}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{2}.tcon.name = 'FPA cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{2}.tcon.convec = repmat([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{2}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{3}.tcon.name = 'FPA target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{3}.tcon.convec = repmat([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{3}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{4}.tcon.name = 'FPA target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{4}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{4}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{5}.tcon.name = 'FPA feedback correct';
            %TDD
            jobs{1}.stats{1}.con.consess{5}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{5}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{6}.tcon.name = 'FPA feedback incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{6}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{6}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{7}.tcon.name = 'NSWP cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{7}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{7}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{8}.tcon.name = 'NSWP cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{8}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{8}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{9}.tcon.name = 'NSWP target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{9}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{9}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{10}.tcon.name = 'NSWP target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{10}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{10}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{11}.tcon.name = 'NSWP feedback correct';
            %TDD
            jobs{1}.stats{1}.con.consess{11}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{11}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{12}.tcon.name = 'NSWP feedback incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{12}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{12}.tcon.sessrep = 'none';
            
            
        elseif strcmp(Task_A{s0},'N')==1 % if NSWP is fist task, append files such that NSWP files are first
            
            jobs{1}.stats{1}.con.consess{1}.tcon.name = 'FPA cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{1}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{1}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{2}.tcon.name = 'FPA cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{2}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{2}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{3}.tcon.name = 'FPA target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{3}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{3}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{4}.tcon.name = 'FPA target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{4}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{4}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{5}.tcon.name = 'FPA feedback correct';
            %TDD
            jobs{1}.stats{1}.con.consess{5}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{5}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{6}.tcon.name = 'FPA feedback incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{6}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{6}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{7}.tcon.name = 'NSWP cue correct';
            %TDD
            jobs{1}.stats{1}.con.consess{7}.tcon.convec = repmat([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{7}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{8}.tcon.name = 'NSWP cue incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{8}.tcon.convec = repmat([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{8}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{9}.tcon.name = 'NSWP target correct';
            %TDD
            jobs{1}.stats{1}.con.consess{9}.tcon.convec = repmat([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{9}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{10}.tcon.name = 'NSWP target incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{10}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{10}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{11}.tcon.name = 'NSWP feedback correct';
            %TDD
            jobs{1}.stats{1}.con.consess{11}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{11}.tcon.sessrep = 'none';
            
            jobs{1}.stats{1}.con.consess{12}.tcon.name = 'NSWP feedback incorrect';
            %TDD
            jobs{1}.stats{1}.con.consess{12}.tcon.convec = repmat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1 n_sess]);
            jobs{1}.stats{1}.con.consess{12}.tcon.sessrep = 'none';
            
            
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