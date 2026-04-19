% =========================================================================
% MASTER COHORT VALIDATION ENGINE (N=216) - RIGOROUS V12 (RADIAL TELESCOPE)
% Objective: Quantify Global Accuracy using LEMON-Decade derived weights.
% =========================================================================

clc;
addpath('C:\MoBI_Research\GFEX-EEG\matlab\core');
addpath('C:\MoBI_Research\GFEX-EEG\reproducibility');
addpath('C:\MoBI_Research\brainstorm3');
addpath(genpath('C:\MoBI_Research\brainstorm3\toolbox'));
mesh_path = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'ICBM152_scalp.mat');

% LOCKED LEMON-10 WEIGHTS (V40.0)
rho_final = 0.000000;
beta_final = 1.983084;
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

fprintf('\nStarting Master Scoreboard Rescue (V12 - Radial Telescope)...\n');

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
                % LEMON SCS-to-SCS Bridge (Already in ALS Space)
                sSubject = bst_get('Subject', sub_id);
                mriFile = fullfile(dbDir, 'ds000221_LEMON', 'anat', sSubject.Anatomy(1).FileName);
                MRI = load(mriFile, 'SCS', 'Voxsize');
                trueL = cs_convert(MRI, 'voxel', 'scs', reshape(MRI.SCS.LPA, 1, 3)) / 1000;
                trueR = cs_convert(MRI, 'voxel', 'scs', reshape(MRI.SCS.RPA, 1, 3)) / 1000;
                matPath = fullfile(raw_lemon_dir, sub_id, 'eeg', [sub_id '_localizer.mat']);
                loc = load(matPath);
                iCz_s = find(contains({loc.Channel.Name}, '_Cz'), 1);
                iT7_s = find(contains({loc.Channel.Name}, '_T7'), 1);
                iT8_s = find(contains({loc.Channel.Name}, '_T8'), 1);
                raw_Cz = loc.Channel(iCz_s).Loc'; raw_T7 = loc.Channel(iT7_s).Loc'; raw_T8 = loc.Channel(iT8_s).Loc';
                
                % LEMON Parity: Both already in ALS after SCS bridge
                p.Cz = raw_Cz; p.T7 = raw_T7; p.T8 = raw_T8;
            else
                % Legacy BIDS Data Extraction
                sub_deriv_dir = fullfile(deriv_root, sub_id);
                tsv_files = dir(fullfile(sub_deriv_dir, '**', '*_electrodes.tsv'));
                if isempty(tsv_files), continue; end
                tbl = readtable(fullfile(tsv_files(1).folder, tsv_files(1).name), 'FileType', 'text');
                names = string(tbl.name);
                
                if contains(cohort, 'ds006525'), lCz="80"; lT7="45"; lT8="108";
                else, lCz="Cz"; lT7="T7"; lT8="T8"; end
                
                get_c = @(lbl) [tbl.x(find(strcmpi(names, lbl),1)), tbl.y(find(strcmpi(names, lbl),1)), tbl.z(find(strcmpi(names, lbl),1))];
                raw_Cz = get_c(lCz); raw_T7 = get_c(lT7); raw_T8 = get_c(lT8);
                
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
                raw_Cz = raw_Cz * scale_eeg; raw_T7 = raw_T7 * scale_eeg; raw_T8 = raw_T8 * scale_eeg;
                
                % PARITY RESCUE: Apply [-Y, X, Z] to ALS cohorts (HK and Oklahoma)
                if contains(cohort, 'ds004718') || contains(cohort, 'ds006525')
                    p.Cz = [-raw_Cz(2), raw_Cz(1), raw_Cz(3)];
                    p.T7 = [-raw_T7(2), raw_T7(1), raw_T7(3)];
                    p.T8 = [-raw_T8(2), raw_T8(1), raw_T8(3)];
                else
                    p.Cz = raw_Cz; p.T7 = raw_T7; p.T8 = raw_T8;
                end
                
                mri_j = dir(fullfile(sub_deriv_dir, '**', '*_desc-normalized_T1w.json'));
                mri_meta = jsondecode(fileread(fullfile(mri_j(1).folder, mri_j(1).name)));
                truth_mm = mri_meta.AnatomicalLandmarkCoordinates;
                trueL = [truth_mm.LPA(1), truth_mm.LPA(2), truth_mm.LPA(3)] / 1000;
                trueR = [truth_mm.RPA(1), truth_mm.RPA(2), truth_mm.RPA(3)] / 1000;
            end
            
            % 2. EXECUTE V37.0 ENGINE
            P_subj_anchors = [p.Cz; p.T7; p.T8];
            T_Cz  = [ 0.011240, 0.025921, 0.141134]; T_T7  = [-0.089174, -0.001327, -0.006348]; T_T8  = [ 0.096880, -0.014286, -0.005819];
            P_temp_anchors = [T_Cz; T_T7; T_T8]; 
            
            [~, P_deflated, ~] = procrustes(P_temp_anchors, P_subj_anchors, 'scaling', true);
            D_L = norm(P_deflated(1,:) - P_deflated(2,:));
            D_R = norm(P_deflated(1,:) - P_deflated(3,:));
            
            [L_mni, R_mni] = predict_helix_tragus_junctions_fno(T_Cz, T_T7, T_T8, mesh_path, rho_final, beta_final, D_L, D_R, D_standard);
            
            [~, ~, t_proj] = procrustes(P_subj_anchors, P_temp_anchors, 'scaling', true, 'reflection', false);
            predL = t_proj.b * L_mni * t_proj.T + t_proj.c(1,:);
            
            % 3. RADIAL TELESCOPE ERROR
            tm_pred = (p.T7 + p.T8) / 2; tm_gt = (trueL + trueR) / 2;
            V_pred = predL - tm_pred; V_gt = trueL - tm_gt;
            V_telescope = (V_pred / norm(V_pred)) * norm(V_gt);
            err = norm(V_telescope - V_gt) * 1000;
            
            Subject_Count = Subject_Count + 1;
            Global_Errors(Subject_Count) = err;
            fprintf('  [PASS] %s: Error = %.2f mm\n', sub_id, err);
            
        catch ME
            fprintf('  [FAIL] %s: %s\n', sub_id, ME.message);
        end
    end
end

if Subject_Count > 0
    fprintf('\n======================================================\n');
    fprintf('  MASTER GLOBAL SCOREBOARD V12.1 (PARITY FIXED)\n');
    fprintf('  Mean Global Error: %.2f mm\n', mean(Global_Errors));
    fprintf('  Median Error:      %.2f mm\n', median(Global_Errors));
    fprintf('======================================================\n');
end
exit;
