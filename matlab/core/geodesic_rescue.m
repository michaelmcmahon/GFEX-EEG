function [pLHJ, pRHJ, info] = geodesic_rescue(varargin)
%GEODESIC_RESCUE  Predict left/right helix-tragus junction from three scalp anchors.
%
%   [pLHJ, pRHJ, info] = GEODESIC_RESCUE(Cz, T7, T8) runs the deterministic
%   geodesic walk (V2.2) on the ICBM152 scalp manifold and returns the
%   predicted helix-tragus junction coordinates (LHJ, RHJ) in
%   right-anterior-superior (RAS) metres. INFO contains diagnostic fields
%   from the walk (path arc lengths, scaling ratios, anchors used).
%
%   USAGE
%       % Default LEMON-tuned weights
%       [pLHJ, pRHJ] = geodesic_rescue(Cz, T7, T8);
%
%       % Cohort preset (loads preset rho / beta / D_standard)
%       [pLHJ, pRHJ] = geodesic_rescue(Cz, T7, T8, 'cohort', 'LEMON_Polhemus_Adult');
%
%       % Opt-in MLP residual correction (Tier 1.5)
%       [pLHJ, pRHJ] = geodesic_rescue(Cz, T7, T8, ...
%                                     'cohort', 'LEMON_Polhemus_Adult', ...
%                                     'mlp_correction', true);
%
%   INPUTS expect metres in RAS by default. Engine is strict meter-space —
%   no internal unit scaling. Set 'parity', 'ALS' for ALS-frame inputs.
%
%   NOTES
%       Production weights (LEMON-tuned 2026-04-19): rho=0.248383,
%       beta=0.235926, D_standard=0.1388. Bit-identical parity with the
%       Python implementation (28 of 28 cross-platform tests at 0.0000 mm
%       wrapper-parity).
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Core extrapolation engine (MATLAB)
% ------------------------------------------------------------------------------
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
%   Issues     : https://github.com/michaelmcmahon/GFEX-EEG/issues
%
%   CITATION (please cite both)
%     [Software] McMahon, M., Schukat, M., & Barrett, E. (2026).
%                GFEX-EEG Toolbox [Software].
%                Zenodo. https://doi.org/10.5281/zenodo.20580899
%     [Paper]    McMahon, M., Schukat, M., & Barrett, E. (Submitted).
%                GFEX-EEG: Geodesic recovery of anatomical fiducials for
%                MRI-free EEG source imaging.
%     See also CITATION.cff in the repository root (machine-readable).
%
%   VERSION
%     For the current package version, call gfex_version() (single source
%     of truth in matlab/core/gfex_version.m).
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
%     See LICENSE in the repository root.
% ==============================================================================

    % 1. Parse Inputs (Flexible Signature)
    if istable(varargin{1}) || isstruct(varargin{1})
        coords = varargin{1};
        [Cz_m, T7_m, T8_m] = extract_standard_anchors(coords);
        offset = 2;
    else
        Cz_m = varargin{1}; T7_m = varargin{2}; T8_m = varargin{3};
        offset = 4;
    end
    
    % Optional Params (hard defaults mirror weight_zoo.json 'default' alias)
    rho = 0.248383; beta = 0.235926; mesh_path = ''; parity = 'RAS'; D_standard = 0.1388;
    cohort = '';
    mlp_correction = false;   % opt-in MLP residual correction (requires cohort)
    rho_set = false; beta_set = false; Dstd_set = false;
    for i = offset:2:nargin
        if strcmpi(varargin{i}, 'rho'),   rho = varargin{i+1}; rho_set = true; end
        if strcmpi(varargin{i}, 'beta'),  beta = varargin{i+1}; beta_set = true; end
        if strcmpi(varargin{i}, 'mesh'),  mesh_path = varargin{i+1}; end
        if strcmpi(varargin{i}, 'parity'), parity = varargin{i+1}; end
        if strcmpi(varargin{i}, 'D_standard'), D_standard = varargin{i+1}; Dstd_set = true; end
        if strcmpi(varargin{i}, 'cohort'), cohort = varargin{i+1}; end
        if strcmpi(varargin{i}, 'mlp_correction'), mlp_correction = varargin{i+1}; end
    end

    % Cohort preset lookup (explicit rho/beta/D_standard kwargs win)
    if ~isempty(cohort)
        preset = load_cohort_preset(cohort);
        if ~rho_set,  rho  = preset.rho;        end
        if ~beta_set, beta = preset.beta;       end
        if ~Dstd_set, D_standard = preset.D_standard; end
    end
    
    if strcmpi(parity, 'ALS') || strcmpi(parity, 'BIDS-Brainstorm')
        % Apply [-Y, X, Z] swap to restore native RAS
        Cz_m = [-Cz_m(2), Cz_m(1), Cz_m(3)];
        T7_m = [-T7_m(2), T7_m(1), T7_m(3)];
        T8_m = [-T8_m(2), T8_m(1), T8_m(3)];
    end
    
    if isempty(mesh_path)
        base_path = fileparts(mfilename('fullpath'));
        mesh_path = fullfile(base_path, '..', '..', 'data', 'ICBM152_scalp.mat');
    end
    
    % 2. Execute Geodesic Engine (MNI Space Prediction)
    % Template Anchors (Golden Kite)
    T_Cz = [0.011240, 0.025921, 0.141134]; 
    T_T7 = [-0.089174, -0.001327, -0.006348]; 
    T_T8 = [0.096880, -0.014286, -0.005819];
    
    D_L_sub = norm(Cz_m - T7_m);
    D_R_sub = norm(Cz_m - T8_m);
    
    % Call Core Engine (predict_helix_tragus_junctions_fno)
    [Lp_mni, Rp_mni, iCz, iT7, iT8] = predict_helix_tragus_junctions_fno(T_Cz, T_T7, T_T8, mesh_path, rho, beta, D_L_sub, D_R_sub, D_standard);
    
    % 3. Temporal Pivot Mapping (Subject Space)
    mesh_data = load(mesh_path, 'Vertices');
    V_temp = mesh_data.Vertices;
    
    P_temp = [V_temp(iCz,:); V_temp(iT7,:); V_temp(iT8,:)];
    P_subj = [Cz_m; T7_m; T8_m];

    [~, ~, tr] = procrustes(P_subj, P_temp, 'scaling', true, 'reflection', false);
    
    t_pivot = (T_T7 + T_T8) / 2;
    s_pivot = (T7_m + T8_m) / 2;
    
    pLHJ = tr.b * (Lp_mni(1,:) - t_pivot) * tr.T + s_pivot;
    pRHJ = tr.b * (Rp_mni(1,:) - t_pivot) * tr.T + s_pivot;
    
    info.rho = rho;
    info.beta = beta;
    info.D_standard = D_standard;
    info.cohort = cohort;
    info.procrustes_residual = norm(P_subj - (tr.b * P_temp * tr.T + tr.c(1,:)));

    % 4. OPT-IN MLP RESIDUAL CORRECTION (Tier 1.5)
    % Applied after pivot-lock projection. Loads cohort-specific MLP weights
    % file per weight_zoo.json 'mlp_weights_file' entry; if absent, errors.
    info.mlp_applied = false;
    if mlp_correction
        if isempty(cohort)
            error('GeodesicRescue:MLPCohortRequired', ...
                'mlp_correction=true requires a cohort kwarg for weight lookup.');
        end
        preset = load_cohort_preset(cohort);
        if ~isfield(preset, 'mlp_weights_file') || isempty(preset.mlp_weights_file)
            error('GeodesicRescue:MLPNotAvailable', ...
                'Cohort "%s" has no mlp_weights_file registered in weight_zoo.json.', cohort);
        end
        base_path = fileparts(mfilename('fullpath'));
        mlp_path = fullfile(base_path, '..', '..', 'data', 'mlp', preset.mlp_weights_file);
        [pLHJ, pRHJ, mlp_info] = apply_mlp_correction(Cz_m, T7_m, T8_m, pLHJ, pRHJ, mlp_path);
        info.mlp_applied    = true;
        info.mlp_weights    = preset.mlp_weights_file;
        info.mlp_dL         = mlp_info.dL;
        info.mlp_dR         = mlp_info.dR;
        info.mlp_correction_magnitude_L_mm = mlp_info.dL_magnitude_mm;
        info.mlp_correction_magnitude_R_mm = mlp_info.dR_magnitude_mm;
    end
end

function [Cz, T7, T8] = extract_standard_anchors(coords)
    if istable(coords)
        labels = string(coords.name);
        x = coords.x; y = coords.y; z = coords.z;
    else
        labels = string(coords.label);
        x = coords.pos(:,1); y = coords.pos(:,2); z = coords.pos(:,3);
    end
    find_pt = @(targets) find_by_label(labels, targets, x, y, z);
    Cz = find_pt({'Cz', 'CZ', '80', '18', '36', 'E36', 'EEG034'});
    T7 = find_pt({'T7', 'T3', '45', '13', '44', '46', 'E45', 'E46', 'EEG030'});
    T8 = find_pt({'T8', 'T4', '108', '14', '77', '102', 'E108', 'E102', 'EEG038'});
    if isempty(Cz) || isempty(T7) || isempty(T8)
        error('GeodesicRescue:MissingAnchors', 'Could not find Cz, T7, or T8 labels.');
    end
end

function pos = find_by_label(labels, targets, x, y, z)
    pos = [];
    for i = 1:length(labels)
        lbl = labels(i);
        stripped = regexprep(lbl, '^(E|EEG)', '', 'ignorecase');
        if any(strcmpi(lbl, targets)) || any(strcmpi(stripped, targets))
            xi = x(i); if iscell(xi), xi = xi{1}; end
            yi = y(i); if iscell(yi), yi = yi{1}; end
            zi = z(i); if iscell(zi), zi = zi{1}; end
            pos = [double(xi), double(yi), double(zi)];
            return;
        end
    end
end
