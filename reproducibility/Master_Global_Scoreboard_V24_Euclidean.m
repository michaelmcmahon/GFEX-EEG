% ==============================================================================
%   GFEX-EEG TOOLBOX — Master cohort validation scoreboard (V24, MATLAB)
% ------------------------------------------------------------------------------
%   MASTER COHORT VALIDATION ENGINE (N=216) - RIGOROUS V24 (PURE EUCLIDEAN)
%   Objective: Quantify Global Accuracy using un-masked 3D Euclidean Distance.
%
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
%
%   CITATION (please cite both)
%     [Software] McMahon, M., Schukat, M., & Barrett, E. (2026).
%                GFEX-EEG Toolbox [Software].
%                Zenodo. https://doi.org/10.5281/zenodo.20580899
%     [Paper]    McMahon, M., Schukat, M., & Barrett, E. (Submitted).
%                GFEX-EEG: Geodesic recovery of anatomical fiducials for
%                MRI-free EEG source imaging.
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
% ==============================================================================

clc;
addpath('C:\MoBI_Research\GFEX-EEG\matlab\core');
addpath('C:\MoBI_Research\GFEX-EEG\reproducibility');
addpath('C:\MoBI_Research\brainstorm3');
addpath(genpath('C:\MoBI_Research\brainstorm3\toolbox'));
mesh_path = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'ICBM152_scalp.mat');

% CROSS-VALIDATED UNIVERSAL CONSTANTS (NO MASKING)
rho_final = 0.248383;
beta_final = 0.235926;
D_standard = 0.1388;

base_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp\Data_Clean';
raw_lemon_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp\Data_Raw\ds000221_LEMON';

% Target Cohorts
cohorts = {'ds000221_LEMON', 'ds004718', 'ds005795', 'ds005811-ds005810', 'ds006525', 'ds007216', 'ds007353'};

Results = struct();
Subject_Count = 0;
Global_Errors = [];

brainstorm server;
dbDir = bst_get('BrainstormDbDir');

fprintf('\nStarting Master Scoreboard Rescue (V24 - UN-MASKED EUCLIDEAN)...\n');

for c = 1:length(cohorts)
    cohort = cohorts{c};
    fprintf('\n--- Processing Cohort: %s ---\n', cohort);

    if strcmpi(cohort, 'ds000221_LEMON')
        sub_ids = {'sub-010060', 'sub-010061', 'sub-010062', 'sub-010063', 'sub-010064', ...
                   'sub-010065', 'sub-010066', 'sub-010067', 'sub-010068', 'sub-010069'};
        iProtocol = bst_get('Protocol', 'ds000221_LEMON');
        gui_brainstorm('SetCurrentProtocol', iProtocol);
    else
        deriv_root = fullfile(base_dir, cohort, 'derivatives', 'intensity_normalization');
        if ~exist(deriv_root, 'dir'), continue; end
        items = dir(deriv_root);
        sub_ids = {items([items.isdir] & startsWith({items.name}, 'sub-')).name};
    end

    for s = 1:length(sub_ids)
        sub_id = sub_ids{s};
        try
            if strcmpi(cohort, 'ds000221_LEMON')
                sSubject = bst_get('Subject', sub_id);
                mriFile = fullfile(dbDir, 'ds000221_LEMON', 'anat', sSubject.Anatomy(1).FileName);
                MRI = load(mriFile, 'SCS', 'Voxsize');
                
                % INTELLIGENT EXTRACTION: Using cs_convert to bridge any Z-offset anomalies
                trueL = cs_convert(MRI, 'voxel', 'scs', reshape(MRI.SCS.LPA, 1, 3)) / 1000;
                trueR = cs_convert(MRI, 'voxel', 'scs', reshape(MRI.SCS.RPA, 1, 3)) / 1000;
                
                matPath = fullfile(raw_lemon_dir, sub_id, 'eeg', [sub_id '_localizer.mat']);
                loc = load(matPath);
                iCz_s = find(contains({loc.Channel.Name}, '_Cz'), 1);
                iT7_s = find(contains({loc.Channel.Name}, '_T7'), 1);
                iT8_s = find(contains({loc.Channel.Name}, '_T8'), 1);
                p.Cz = loc.Channel(iCz_s).Loc'; p.T7 = loc.Channel(iT7_s).Loc'; p.T8 = loc.Channel(iT8_s).Loc';
            else
                % Legacy BIDS Data Extraction (with RAS Parity Restoration)
                sub_deriv_dir = fullfile(deriv_root, sub_id);
                tsv_files = dir(fullfile(sub_deriv_dir, '**', '*_electrodes.tsv'));
                if isempty(tsv_files), continue; end
                tbl = readtable(fullfile(tsv_files(1).folder, tsv_files(1).name), 'FileType', 'text');
                names = string(tbl.name);
                if contains(cohort, 'ds006525'), lCz="80"; lT7="45"; lT8="108";
                else, lCz="Cz"; lT7="T7"; lT8="T8"; end
                get_c = @(lbl) [tbl.x(find(strcmpi(names, lbl),1)), tbl.y(find(strcmpi(names, lbl),1)), tbl.z(find(strcmpi(names, lbl),1))];
                raw_Cz = get_c(lCz); raw_T7 = get_c(lT7); raw_T8 = get_c(lT8);
                
                % Handle scaling (micrometer trap detection)
                scale_eeg = 1.0;
                cs_files = dir(fullfile(sub_deriv_dir, '**', '*_coordsystem.json'));
                if ~isempty(cs_files)
                    cs = jsondecode(fileread(fullfile(cs_files(1).folder, cs_files(1).name)));
                    if isfield(cs, 'EEGCoordinateUnits')
                        if strcmpi(cs.EEGCoordinateUnits, 'cm'), scale_eeg = 0.01;
                        elseif strcmpi(cs.EEGCoordinateUnits, 'mm'), scale_eeg = 0.001;
                        end
                    end
                end
                p.Cz = raw_Cz * scale_eeg; p.T7 = raw_T7 * scale_eeg; p.T8 = raw_T8 * scale_eeg;

                % Apply ALS-to-RAS parity rescue
                if contains(cohort, 'ds004718') || contains(cohort, 'ds006525') || contains(cohort, 'ds007353')
                    p.Cz = [-p.Cz(2), p.Cz(1), p.Cz(3)];
                    p.T7 = [-p.T7(2), p.T7(1), p.T7(3)];
                    p.T8 = [-p.T8(2), p.T8(1), p.T8(3)];
                end
                
                mri_j = dir(fullfile(sub_deriv_dir, '**', '*_desc-normalized_T1w.json'));
                mri_meta = jsondecode(fileread(fullfile(mri_j(1).folder, mri_j(1).name)));
                truth_mm = mri_meta.AnatomicalLandmarkCoordinates;
                trueL = [truth_mm.LPA(1), truth_mm.LPA(2), truth_mm.LPA(3)] / 1000;
            end

            % 2. EXECUTE PRODUCTION ENGINE (V22.2 Locked)
            [predL, predR, info] = geodesic_rescue(p.Cz, p.T7, p.T8, 'rho', rho_final, 'beta', beta_final, 'mesh', mesh_path);

            % 3. UN-MASKED PURE EUCLIDEAN ERROR
            err = norm(predL - trueL) * 1000;

            Subject_Count = Subject_Count + 1;
            Global_Errors(Subject_Count) = err;
            fprintf('  [PASS] %s: TRE = %.2f mm\n', sub_id, err);

        catch ME
            fprintf('  [FAIL] %s: %s\n', sub_id, ME.message);
        end
    end
end

if Subject_Count > 0
    fprintf('\n======================================================\n');
    fprintf('  MASTER GLOBAL SCOREBOARD V24.0 (UN-MASKED)\n');
    fprintf('  Mean Global TRE:   %.2f mm\n', mean(Global_Errors));
    fprintf('  Median TRE:        %.2f mm\n', median(Global_Errors));
    fprintf('======================================================\n');
end
exit;
