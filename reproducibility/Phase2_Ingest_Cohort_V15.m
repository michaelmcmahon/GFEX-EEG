function Phase2_Ingest_Cohort_V15(cohort_id)
%PHASE2_INGEST_COHORT_V15  Universal BIDS cohort ingestion with fiducial capture.
%
%   Phase 2 universal ingestion (V15.1, fiducial-enabled). Ingests any
%   cohort from the canonical BIDS layout and captures EEG fiducials
%   ready for the V19 origin-bridging Procrustes step. Outputs are
%   written to ``Data_Processed/<cohort_id>/`` and the matching
%   Brainstorm protocol under ``brainstorm_db/<cohort_id>/anat/``.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Phase 2 universal BIDS ingestion (MATLAB, V15.1)
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

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    data_dir = fullfile(root_dir, 'Data_Clean', cohort_id);
    proc_dir = fullfile(root_dir, 'Data_Processed', cohort_id);
    bst_dir  = fullfile('C:\MoBI_Research\brainstorm_db', cohort_id, 'anat');
    
    if ~exist(proc_dir, 'dir'), mkdir(proc_dir); end
    
    scripts_dir = fileparts(mfilename('fullpath'));
    addpath(scripts_dir);
    
    deriv_root = fullfile(data_dir, 'derivatives', 'intensity_normalization');
    if ~exist(deriv_root, 'dir'), deriv_root = data_dir; end
    
    s_dirs = dir(fullfile(deriv_root, 'sub-*'));
    s_dirs = s_dirs([s_dirs.isdir]);
    
    Master_Scalp_Proxy = struct();
    Ground_Truth_MRI = struct();
    
    fprintf('Starting V15.1 Ingestion for %s (%d subjects)...\n', cohort_id, length(s_dirs));
    
    for i = 1:length(s_dirs)
        sub_id_raw = s_dirs(i).name;
        sub_id = strrep(sub_id_raw, '-', '_');
        fprintf('  Processing %s... ', sub_id_raw);
        
        try
            % 1. EEG SENSORS & FIDUCIALS
            tsv_files = dir(fullfile(deriv_root, sub_id_raw, '**', '*_electrodes.tsv'));
            valid_tsv = [];
            for j = 1:length(tsv_files)
                if ~startsWith(tsv_files(j).name, '._'), valid_tsv = [valid_tsv, j]; end
            end
            if isempty(valid_tsv), error('TSV missing'); end
            tsv_path = fullfile(tsv_files(valid_tsv(1)).folder, tsv_files(valid_tsv(1)).name);
            
            % Sanitize sensors
            Proxy = universal_bids_sanitizer_v15(cohort_id, tsv_path);
            
            % Extract EEG Fiducials from coordsystem.json
            json_files = dir(fullfile(deriv_root, sub_id_raw, '**', '*_coordsystem.json'));
            valid_json = [];
            for j = 1:length(json_files)
                if ~startsWith(json_files(j).name, '._'), valid_json = [valid_json, j]; end
            end
            
            if ~isempty(valid_json)
                json_path = fullfile(json_files(valid_json(1)).folder, json_files(valid_json(1)).name);
                cs = jsondecode(fileread(json_path));
                
                % Standardize units (trap check)
                scale = 1.0;
                if isfield(cs, 'EEGCoordinateUnits')
                    if strcmpi(cs.EEGCoordinateUnits, 'mm'), scale = 0.001;
                    elseif strcmpi(cs.EEGCoordinateUnits, 'cm'), scale = 0.01;
                    elseif strcmpi(cs.EEGCoordinateUnits, 'm'), scale = 1.0;
                    end
                end
                
                % Capture LPA, RPA, NAS
                if isfield(cs, 'AnatomicalLandmarkCoordinates')
                    f = cs.AnatomicalLandmarkCoordinates;
                    nas = []; if isfield(f, 'NAS'), nas = f.NAS; elseif isfield(f, 'Nasion'), nas = f.Nasion; end
                    lpa = []; if isfield(f, 'LPA'), lpa = f.LPA; end
                    rpa = []; if isfield(f, 'RPA'), rpa = f.RPA; end
                    
                    if ~isempty(nas) && ~isempty(lpa) && ~isempty(rpa)
                        % Apply V15 Scaling/Rotation to Fiducials if needed (e.g. HK, OK, WH)
                        if strcmp(cohort_id, 'ds004718')
                            trans = @(v) [-v(2), v(1), v(3)] * 0.001;
                            Proxy.LPA = trans(lpa); Proxy.RPA = trans(rpa); Proxy.NAS = trans(nas);
                        elseif strcmp(cohort_id, 'ds006525') || strcmp(cohort_id, 'ds000117-ds002718')
                            % ALS to RAS for OK and WH
                            trans = @(v) [-v(2), v(1), v(3)] * scale;
                            Proxy.LPA = trans(lpa); Proxy.RPA = trans(rpa); Proxy.NAS = trans(nas);
                        else
                            Proxy.LPA = lpa(:)' * scale; Proxy.RPA = rpa(:)' * scale; Proxy.NAS = nas(:)' * scale;
                        end
                    end
                end
            end
            
            Master_Scalp_Proxy.(sub_id) = Proxy;
            
            % 2. MRI GROUND TRUTH
            bst_sub_dir = fullfile(bst_dir, [cohort_id, '_', sub_id_raw]);
            if ~exist(bst_sub_dir, 'dir') && strcmp(cohort_id, 'ds007216')
                bst_sub_dir = fullfile(bst_dir, ['ds00721601_', sub_id_raw]);
            end
            if ~exist(bst_sub_dir, 'dir'), bst_sub_dir = fullfile(bst_dir, sub_id_raw); end
            
            mri_files = dir(fullfile(bst_sub_dir, 'subjectimage_*.mat'));
            if isempty(mri_files), error('Brainstorm MRI missing'); end
            mri_path = fullfile(mri_files(1).folder, mri_files(1).name);
            Ground_Truth_MRI.(sub_id).mri_path = mri_path;
            
            fprintf('OK\n');
        catch ME
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    
    save(fullfile(proc_dir, 'Master_Scalp_Proxy.mat'), 'Master_Scalp_Proxy');
    save(fullfile(proc_dir, 'Ground_Truth_MRI.mat'), 'Ground_Truth_MRI');
end
