% Pfad zur SPM-Direktory und den Daten
basePath = '/Volumes/CJ_NetzTran/fMRI';
conPath_IR = fullfile(basePath, '1st_Level/Memory_Immediate_Recall');
conPath_DR = fullfile(basePath, '1st_Level/Memory_Delayed_Recall');
% resultsDir = fullfile(basePath, '2nd_Level/2x2Anova_TIME_TASK_APM_LGT_cov');
resultsDir = fullfile(basePath, '2nd_Level/2x2Anova_TIME_TASK_APM_LGT_cov_29');
addpath('/Applications/MATLAB_R2021b.app/toolbox/spm12')

% Sicherstellen, dass das Ergebnisverzeichnis existiert
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Liste der eingeschlossenen Proband*innen
subjects = {'VP02','VP03',...};
% Initialisieren der matlabbatch-Variable
matlabbatch = {};

% Spezifizieren des Ausgabeordners für das faktorenielle Design
matlabbatch{1}.spm.stats.factorial_design.dir = {resultsDir};

% Konfigurieren des faktoriellen Designs
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Timepoints';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % Abhängige Messungen
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1; % Ungleiche Varianz
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0; % Keine Skalierung des Gesamtmittels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0; % Keine Kovariate

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'StimulusType';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1; % Abhängige Messungen
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1; % Ungleiche Varianz
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0; % Keine Skalierung des Gesamtmittels
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0; % Keine Kovariate

% Initialisieren des iCell-Zählers und des Scans-Arrays
iCell = 0;
allScans = {};
numSubjects = length(subjects);

% Definieren der Zellen und Zuweisen der Kontrastdateien
for t = 1:2  % Zeitpunkte: 1-Immediate, 2-Delayed
    for st = 1:2  % Stimulus Typ: 1-FPA, 2-NSWP
        scansForCell = {};
        for s = 1:numSubjects
            if t == 1 && st == 1
                conFile = fullfile(conPath_IR, subjects{s}, 'con_0001.nii'); % Immediate Recall, FPA
            elseif t == 1 && st == 2
                conFile = fullfile(conPath_IR, subjects{s}, 'con_0004.nii'); % Immediate Recall, NSWP
            elseif t == 2 && st == 1
                conFile = fullfile(conPath_DR, subjects{s}, 'con_0001.nii'); % Delayed Recall, FPA
            elseif t == 2 && st == 2
                conFile = fullfile(conPath_DR, subjects{s}, 'con_0003.nii'); % Delayed Recall, NSWP
            end
            
            % Überprüfen, ob die Datei existiert
            if exist(conFile, 'file')
                scansForCell{end+1} = conFile;
            else
                warning('Datei %s existiert nicht. Überspringe...', conFile);
            end
        end
        
        if ~isempty(scansForCell)
            iCell = iCell + 1;
            matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(iCell).levels = [t, st];
            matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(iCell).scans = scansForCell';
            allScans = [allScans, scansForCell];
        end
    end
end

% Sicherstellen, dass mindestens ein Scan gefunden wurde
if isempty(allScans)
    error('Keine gültigen Kontrastdateien gefunden. Bitte überprüfen Sie die Pfade und Dateinamen.');
end
% Hinzufügen der Probanden-Kovariaten
subjectIndex = 0;
for s = 1:numSubjects
    covariateVector = zeros(length(allScans), 1);
    for i = 1:length(allScans)
        if contains(allScans{i}, subjects{s})
            covariateVector(i) = 1;
        end
    end
    if sum(covariateVector) > 0
        subjectIndex = subjectIndex + 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex).c = covariateVector;
        matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex).cname = ['Subject_' subjects{s}];
        matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex).iCFI = 1; % Interaktionen mit Faktoren
        matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex).iCC = 1; % Zentrierung: Gesamtmittel
    end
end

% Laden und Wiederholen der Kovariaten
APM = load('/Volumes/CJ_NetzTran/IQ:MQ/APM29.txt');
LGT = load('/Volumes/CJ_NetzTran/IQ:MQ/LGT29.txt');

APM_cov = repmat(APM, 4, 1); % Wiederhole APM für alle Bedingungen
LGT_cov = repmat(LGT, 4, 1); % Wiederhole LGT für alle Bedingungen

% Definieren der globalen Kovariaten für APM und LGT
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 1).c = APM_cov;
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 1).cname = 'APM';
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 1).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 2).c = LGT_cov;
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 2).cname = 'LGT';
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(subjectIndex + 2).iCC = 1;


% Maskierung und Global Normalisierung
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Modellschätzung
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(resultsDir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Job ausführen
spm('defaults', 'FMRI');
spm_jobman('initcfg');
try
    spm_jobman('run', matlabbatch);
catch ME
    disp('Ein Fehler ist beim Ausführen des Batches aufgetreten:');
    disp(ME.message);
end
