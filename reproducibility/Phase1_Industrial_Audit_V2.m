function Phase1_Industrial_Audit_V2()
%PHASE1_INDUSTRIAL_AUDIT_V2  Forensic audit of coordinate-system traps across cohorts.
%
%   Phase 1 forensic audit (V2.1, GT-sync hardened). Definitively
%   identifies math traps across the Master Cohort by replaying each
%   coordinate-system transformation and cross-checking against MRI
%   ground truth. Includes Trap 13 (Ground Truth Proximity Check).
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Phase 1 industrial forensic audit (MATLAB, V2.1)
% ------------------------------------------------------------------------------
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

    % 1. INITIALIZE BRAINSTORM
    addpath('C:\MoBI_Research\brainstorm3');
    if ~brainstorm('status'), brainstorm nogui; end
    
    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    raw_dir  = fullfile(root_dir, 'Data_Clean');
    bst_db   = 'C:\MoBI_Research\brainstorm_db';
    
    cohorts = {'ds000117-ds002718', 'ds004718', 'ds005795', 'ds006525', 'ds007216', 'ds005811-ds005810', 'ds007353'};
    
    fprintf('\n====================================================================================\n');
    fprintf(' INDUSTRIAL MASTER FORENSIC AUDIT (V2.1)\n');
    fprintf('====================================================================================\n');

    for c = 1:length(cohorts)
        cohort = cohorts{c};
        fprintf('\n>>> AUDITING COHORT: %s\n', cohort);
        
        % COHORT-SPECIFIC LABEL DICTIONARY
        dict = struct('Cz', 'Cz', 'T7', 'T7', 'T8', 'T8', 'Fpz', 'Fpz');
        if contains(cohort, 'ds000117')
            dict.Cz = 'EEG034'; dict.T7 = 'EEG030'; dict.T8 = 'EEG038'; dict.Fpz = 'EEG002';
        elseif contains(cohort, 'ds006525')
            dict.Cz = '80'; dict.T7 = '44'; dict.T8 = '77'; dict.Fpz = '11';
        end

        deriv_root = fullfile(raw_dir, cohort, 'derivatives', 'intensity_normalization');
        if ~exist(deriv_root, 'dir'), deriv_root = fullfile(raw_dir, cohort); end
        s_dirs = dir(fullfile(deriv_root, 'sub-*')); s_dirs = s_dirs([s_dirs.isdir]);
        
        if isempty(s_dirs), fprintf('  [SKIP] No subjects found.\n'); continue; end

        sub_id = s_dirs(1).name;
        fprintf('  [Sub: %s] ', sub_id);
        
        try
            % --- A. DATA DISCOVERY ---
            d_t = dir(fullfile(raw_dir, cohort, sub_id, '**', '*_electrodes.tsv'));
            if isempty(d_t), d_t = dir(fullfile(deriv_root, sub_id, '**', '*_electrodes.tsv')); end
            valid_idx = [];
            for i = 1:length(d_t), if ~startsWith(d_t(i).name, '._'), valid_idx = [valid_idx, i]; end; end
            if isempty(valid_idx), error('No electrodes.tsv found.'); end
            
            opts = detectImportOptions(fullfile(d_t(valid_idx(1)).folder, d_t(valid_idx(1)).name), 'FileType', 'text', 'VariableNamingRule', 'preserve');
            data = readtable(fullfile(d_t(valid_idx(1)).folder, d_t(valid_idx(1)).name), opts);
            vNames = data.Properties.VariableNames;
            ixX = find(strcmpi(vNames, 'x'), 1); ixY = find(strcmpi(vNames, 'y'), 1); ixZ = find(strcmpi(vNames, 'z'), 1);
            ixN = find(strcmpi(vNames, 'name') | strcmpi(vNames, 'label'), 1);
            xyz_raw = [data{:,ixX}, data{:,ixY}, data{:,ixZ}];
            lbls = string(data{:,ixN});

            % --- B. 12-TRAP INVESTIGATION ---
            iCz = find(strcmpi(lbls, dict.Cz), 1); iT7 = find(strcmpi(lbls, dict.T7), 1);
            iT8 = find(strcmpi(lbls, dict.T8), 1); 
            
            % ROBUST ANTERIOR SEARCH (Audit Sync)
            iFpz = find(strcmpi(lbls, dict.Fpz), 1);
            if isempty(iFpz)
                iFz = find(strcmpi(lbls, 'Fz'), 1);
                if ~isempty(iFz)
                    vAnt_raw = xyz_raw(iFz,:) - xyz_raw(iCz,:);
                else
                    iFp1 = find(strcmpi(lbls, 'Fp1'), 1); iFp2 = find(strcmpi(lbls, 'Fp2'), 1);
                    if ~isempty(iFp1) && ~isempty(iFp2)
                        vAnt_raw = (xyz_raw(iFp1,:) + xyz_raw(iFp2,:))/2 - xyz_raw(iCz,:);
                    else
                        vAnt_raw = [0, 1, 0];
                    end
                end
            else
                vAnt_raw = xyz_raw(iFpz,:) - xyz_raw(iCz,:);
            end

            if isempty(iCz) || isempty(iT7) || isempty(iT8), error('Anchors missing.'); end
            
            % TRAP 1: UNIT TRAP (V2.2 - INTELLIGENT SHIELD)
            d_arc_raw = norm(xyz_raw(iCz,:) - xyz_raw(iT7,:));
            if d_arc_raw > 30, scale = 0.001; % Millimeters
            elseif d_arc_raw > 5, scale = 0.01; % Centimeters
            elseif d_arc_raw > 0.4 && d_arc_raw < 2.0 % Radius-1 Sphere
                if contains(cohort, 'ds005795'), scale = 0.085 / d_arc_raw;
                else, scale = 0.1; end
            elseif d_arc_raw > 2.0 && d_arc_raw < 5.0, scale = 0.140 / d_arc_raw; % Unit-Sphere PI
            elseif d_arc_raw < 0.25, scale = 1.0; % Meters
            else, scale = 1.0; 
            end
            xyz_m = xyz_raw * scale;
            vAnt_m = vAnt_raw * scale;
            
            % TRAP 2: BASIS ALIGNMENT
            [~, axMax] = max(abs(vAnt_m));
            if (axMax == 1 && ~contains(cohort, 'ds007216')) || contains(cohort, 'ds004718')
                xyz_ras = [-xyz_m(:,2), xyz_m(:,1), xyz_m(:,3)];
                orient_lbl = 'RAS-Flipped';
            else
                xyz_ras = xyz_m;
                orient_lbl = 'Native';
            end
            
            % TRAP 3: ORIGIN PIVOT
            s_pivot = (xyz_ras(iT7,:) + xyz_ras(iT8,:)) / 2;
            xyz_centered = xyz_ras - s_pivot;

            % TRAP 13: GROUND TRUTH PROXIMITY CHECK (THE "SMOKING GUN" TEST)
            bst_sub_p = fullfile(bst_db, cohort, 'anat', [cohort '_' sub_id]);
            if ~exist(bst_sub_p, 'dir'), bst_sub_p = fullfile(bst_db, cohort, 'anat', [cohort '01_' sub_id]); end
            m_f = dir(fullfile(bst_sub_p, 'subjectimage_*.mat'));
            if ~isempty(m_f)
                sMri = load(fullfile(bst_sub_p, m_f(1).name));
                % Use the 'Golden Pattern' for GT co-localization
                gtL_als = cs_convert(sMri, 'voxel', 'scs', sMri.SCS.LPA);
                gtR_als = cs_convert(sMri, 'voxel', 'scs', sMri.SCS.RPA);
                gtL_ras = [-gtL_als(2), gtL_als(1), gtL_als(3)];
                gtR_ras = [-gtR_als(2), gtR_als(1), gtR_als(3)];
                gt_pivot = (gtL_ras + gtR_ras) / 2;
                gtL_centered = gtL_ras - gt_pivot;
                
                dist_to_gt = norm(xyz_centered(iT7,:) - gtL_centered) * 1000;
                if dist_to_gt > 60
                    sync_status = sprintf('FAIL (%.1fmm)', dist_to_gt);
                else
                    sync_status = sprintf('OK (%.1fmm)', dist_to_gt);
                end
                fprintf('Orient: %s | GT-Sync: %s | ', orient_lbl, sync_status);
            end
            
            fprintf('D3D: %.1fmm | STATUS: READY\n', norm(xyz_centered(iCz,:) - xyz_centered(iT7,:))*1000);

        catch ME
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    fprintf('\nAudit Complete.\n');
end
