% statistical analysis for fMRI data using SPM12
%%% Credits:
% % Author:  Charlotte Jeschina
% Email:    charlotte.jeschina@student.uni-luebeck.de
% Date:     2024
% Institute: University of Luebeck, Institute for clinical and experimental Pharmacology and toxicology, AG Marshall
% Project:  NetzTran

% Pfad zur SPM-Direktory und den Daten
basePath = '/Volumes/CJ_NetzTran/fMRI';
conPath_IR = fullfile(basePath, '1st_Level/Memory_Immediate_Recall');
conPath_DR = fullfile(basePath, '1st_Level/Memory_Delayed_Recall');
resultsDir = fullfile(basePath, '2nd_Level/OneSampleTTest_DR_NSWP');
addpath('/Applications/MATLAB_R2021b.app/toolbox/spm12')

% Sicherstellen, dass das Ergebnisverzeichnis existiert
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Liste der eingeschlossenen Proband*innen
subjects = {'VP02',...};


% Initialisieren der matlabbatch-Variable
matlabbatch = {};

% Spezifizieren des Ausgabeordners für das One-Sample T-Test Design
matlabbatch{1}.spm.stats.factorial_design.dir = {resultsDir};

% Konfigurieren des One-Sample T-Test Designs
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {};

for s = 1:length(subjects)
    conFile = fullfile(conPath_DR, subjects{s}, 'con_0003.nii'); % Delayed Recall, NSWP
%     conFile = fullfile(conPath_DR, subjects{s}, 'con_0004.nii'); % target, Delayed Recall, NSWP
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
