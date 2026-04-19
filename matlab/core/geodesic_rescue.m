% /*******************************************************************************
% * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
% * Version: 1.0.0
% * Repository: https://github.com/michaelmcmahon/GFEX-EEG
% * License:  MIT License
% * Authors: Michael McMahon / University of Galway
% * DOI: [If available]
% * Date: 2026
% *
% * [License Text or link to License file]
% *******************************************************************************/

function [pLHJ, pRHJ, info] = geodesic_rescue(varargin)
% GEODESIC_RESCUE (V21.8 - Pure Core)
% Strict Meter-Space Engine for Geodesic Extrapolation.
% Expects Meters, returns Meters. No internal scaling.

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
    rho_set = false; beta_set = false; Dstd_set = false;
    for i = offset:2:nargin
        if strcmpi(varargin{i}, 'rho'),   rho = varargin{i+1}; rho_set = true; end
        if strcmpi(varargin{i}, 'beta'),  beta = varargin{i+1}; beta_set = true; end
        if strcmpi(varargin{i}, 'mesh'),  mesh_path = varargin{i+1}; end
        if strcmpi(varargin{i}, 'parity'), parity = varargin{i+1}; end
        if strcmpi(varargin{i}, 'D_standard'), D_standard = varargin{i+1}; Dstd_set = true; end
        if strcmpi(varargin{i}, 'cohort'), cohort = varargin{i+1}; end
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
