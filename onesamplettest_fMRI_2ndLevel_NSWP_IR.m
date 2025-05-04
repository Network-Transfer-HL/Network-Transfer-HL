% Pfad zur SPM-Direktory und den Daten
basePath = '/Volumes/CJ_NetzTran/fMRI';
conPath_IR = fullfile(basePath, '1st_Level/Memory_Immediate_Recall');
% resultsDir = fullfile(basePath, '2nd_Level/OneSampleTTest_IR_NSWP');
resultsDir = fullfile(basePath, '2nd_Level/29_OneSampleTTest_IR_NSWP');
addpath('/Applications/MATLAB_R2021b.app/toolbox/spm12')

% Sicherstellen, dass das Ergebnisverzeichnis existiert
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Liste der eingeschlossenen Proband*innen
% subjects = {'VP02','VP03','VP06','VP08','VP09','VP11','VP12','VP14','VP15','VP18','VP19','VP20','VP21','VP24','VP25','VP27','VP28','VP30','VP32','VP34','VP38','VP40','VP41','VP42','VP43','VP44','VP45','VP46','VP47','VP49','VP50','VP52','VP53'};
% subjects = {'VP02','VP03','VP06','VP08','VP09','VP11','VP15','VP18','VP19','VP20','VP21','VP24','VP25','VP27','VP28','VP30','VP32','VP34','VP38','VP40','VP41','VP43','VP44','VP45','VP46','VP47','VP49','VP50','VP52','VP53'};
subjects = {'VP03','VP06','VP08','VP09','VP11','VP15','VP16','VP18','VP19','VP20','VP21','VP24','VP25','VP27','VP28','VP31','VP32','VP34','VP38','VP40','VP41','VP43','VP44','VP45','VP46','VP47','VP50','VP52','VP53'};

% Initialisieren der matlabbatch-Variable
matlabbatch = {};

% Spezifizieren des Ausgabeordners für das One-Sample T-Test Design
matlabbatch{1}.spm.stats.factorial_design.dir = {resultsDir};

% Konfigurieren des One-Sample T-Test Designs
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {};

for s = 1:length(subjects)
    conFile = fullfile(conPath_IR, subjects{s}, 'con_0004.nii'); % cue, Immediate Recall, NSWP
%     conFile = fullfile(conPath_IR, subjects{s}, 'con_0005.nii'); % target, Immediate Recall, NSWP
    if exist(conFile, 'file')
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{end+1, 1} = conFile; % Als Spaltenvektor hinzufügen
    else
        warning('Datei %s existiert nicht. Überspringe...', conFile);
    end
end

% Modellschätzung
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(resultsDir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Kontrastmanager
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(resultsDir, 'SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Positive Effect';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

% Job ausführen
spm('defaults', 'FMRI');
spm_jobman('initcfg');
try
    spm_jobman('run', matlabbatch);
catch ME
    disp('Ein Fehler ist beim Ausführen des Batches aufgetreten:');
    disp(ME.message);
end
