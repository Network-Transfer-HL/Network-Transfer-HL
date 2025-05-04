%%% Credits:
% % Author:   Alexandra Sobczak, M.Sc.
% Email:    alexandra.sobczak@uni-luebeck.de
% Date:     2023-03-14 (YYYY/MM/DD)
% Institute:University of Luebeck, IPSY1, Bunzeck Lab
% Project:  NetzTran
% % Co-Author/ Edited by: Charlotte Jeschina
% Email:   charlotte.jeschina@student.uni-luebeck.de
% Date:     2023-10-12 (YYYY/MM/DD)

% Define the base directory where your SPM data are stored
baseDir = '/Volumes/mehr platz/NetzTran/fMRI/1st_Level';

% Define subjects and tasks
subjects = {'VP02',...};
tasks = {'Memory_Immediate_Recall', 'Memory_Delayed_Recall'};

% Start a loop over tasks and subjects
for t = 1
% for t = 1:length(tasks)
    task = tasks{t};
    for i = 1:length(subjects)
        subject = subjects{i};
        subjectDir = fullfile(baseDir, task, subject);
        
        % Load the SPM.mat file
        spm_mat_file = fullfile(subjectDir, 'SPM.mat');
        
        % Load the SPM.mat into SPM variable
        if exist(spm_mat_file, 'file')
            load(spm_mat_file);
        else
            fprintf('SPM.mat not found for %s in %s\n', subject, task);
            continue;
        end

        % Perform the estimation
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat_file};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

        % Run the job
        spm_jobman('run', matlabbatch);
        clear matlabbatch;

        % Define contrasts
        matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_file};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Correct > Incorrect';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Correct';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Incorrect';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.delete = 1;

        % Run the contrast job
        spm_jobman('run', matlabbatch);
        clear matlabbatch;

        fprintf('Processed %s for task %s\n', subject, task);
    end
end

fprintf('All specified models have been estimated and contrasts defined.\n');
