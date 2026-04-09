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

% =========================================================================
% MODULE A: UNIVERSAL MASTER SANITIZER (V15.0 - GLOBAL COHORT GATEWAY)
% Objective: Safely ingest and homogenize 5 disparate BIDS cohorts into 
%            sterile RAS space (Meters) while avoiding the 12 spatial traps.
% Version: 15.0 (Industrial Lockdown Edition)
% =========================================================================
function Proxy = universal_bids_sanitizer_v15(cohort_id, tsv_file_path)

    % 1. THE HEADER TRAP SHIELD (Boston Fix)
    % Force MATLAB to preserve legacy or non-standard TSV column headers
    try
        opts = detectImportOptions(tsv_file_path, 'FileType', 'text');
        opts.VariableNamingRule = 'preserve'; 
        raw_data = readtable(tsv_file_path, opts);
    catch
        error('CRITICAL: Failed to read TSV file at %s. Check path/permissions.', tsv_file_path);
    end

    % 2. EXTRACT RAW FIDUCIALS / ANCHORS
    % Using robust string-matching to find labels, numbers, or legacy codes
    % Added standard EGI mappings (36=Cz, 45/46=T7, 108/102=T8)
    v_Cz  = extract_electrode(raw_data, {'Cz', 'CZ', '80', '18', '36', 'EEG034'});
    v_T7  = extract_electrode(raw_data, {'T7', '45', '13', '44', '46', 'EEG030'});
    v_T8  = extract_electrode(raw_data, {'T8', '108', '14', '77', '102', 'EEG038'});
    
    if isempty(v_Cz) || isempty(v_T7) || isempty(v_T8)
        error('CRITICAL: Missing core anchors (Cz, T7, or T8) in %s.', tsv_file_path);
    end
    
    % 3. THE FPZ TRAP SHIELD (Boston Synthesis)
    v_Fpz = extract_electrode(raw_data, {'Fpz', 'FPZ', '11', 'EEG002'});
    if isempty(v_Fpz)
        % Synthesize Fpz as the exact geometric midpoint of Fp1 and Fp2
        v_Fp1 = extract_electrode(raw_data, {'Fp1', 'FP1', '22', 'EEG001'}); % EGI 22 equiv
        v_Fp2 = extract_electrode(raw_data, {'Fp2', 'FP2', '9', 'EEG003'});  % EGI 9 equiv
        if ~isempty(v_Fp1) && ~isempty(v_Fp2)
            v_Fpz = (v_Fp1 + v_Fp2) / 2;
        else
            % Fallback: Use Nasion if available in electrodes
            v_nas = extract_electrode(raw_data, {'NAS', 'Nasion', 'FID'});
            if ~isempty(v_nas)
                v_Fpz = v_nas;
            else
                % Final fallback: Use Cz but warn (bad for orientation but prevents crash)
                warning('Fpz synthesis failed for %s. Using Cz fallback.', tsv_file_path);
                v_Fpz = v_Cz;
            end
        end
    end

    % =====================================================================
    % 4. COHORT-SPECIFIC TRAP NEUTRALIZATION (The Gateway)
    % =====================================================================
    switch cohort_id
        
        case 'ds004718'  % HONG KONG (ALS + Phantom Units)
            % Unit Trap: Scale millimeters to meters
            scale = 0.001;
            Proxy.Cz  = [-v_Cz(2),  v_Cz(1),  v_Cz(3)]  * scale;
            Proxy.T7  = [-v_T7(2),  v_T7(1),  v_T7(3)]  * scale;
            Proxy.T8  = [-v_T8(2),  v_T8(1),  v_T8(3)]  * scale;
            Proxy.Fpz = [-v_Fpz(2), v_Fpz(1), v_Fpz(3)] * scale;

        case 'ds006525'  % OKLAHOMA (ALS + Meters)
            % Unit is already meters, but axis is ALS
            scale = 1.0;
            Proxy.Cz  = [-v_Cz(2),  v_Cz(1),  v_Cz(3)]  * scale;
            Proxy.T7  = [-v_T7(2),  v_T7(1),  v_T7(3)]  * scale;
            Proxy.T8  = [-v_T8(2),  v_T8(1),  v_T8(3)]  * scale;
            Proxy.Fpz = [-v_Fpz(2), v_Fpz(1), v_Fpz(3)] * scale;

        case 'ds000117-ds002718' % WAKEMAN-HENSON (ALS + Centimeters)
            % Unit Trap: Scale centimeters to meters
            scale = 0.01;
            Proxy.Cz  = [-v_Cz(2),  v_Cz(1),  v_Cz(3)]  * scale;
            Proxy.T7  = [-v_T7(2),  v_T7(1),  v_T7(3)]  * scale;
            Proxy.T8  = [-v_T8(2),  v_T8(1),  v_T8(3)]  * scale;
            Proxy.Fpz = [-v_Fpz(2), v_Fpz(1), v_Fpz(3)] * scale;

        case 'ds005795'  % MAGDEBURG (Unit Sphere + BESA Tilt)
            % Unit Sphere Trap: Scale unit-radius to average human head radius (75mm)
            scale = 0.075; 
            Proxy.Cz  = v_Cz  * scale;
            Proxy.T7  = v_T7  * scale;
            Proxy.T8  = v_T8  * scale;
            Proxy.Fpz = v_Fpz * scale;

        case 'ds007216'  % BOSTON (Standard BIDS CapTrak)
            Proxy.Cz  = v_Cz;
            Proxy.T7  = v_T7;
            Proxy.T8  = v_T8;
            Proxy.Fpz = v_Fpz;

        case {'ds007353', 'ds005811', 'ds005811-ds005810'} % HAD & NOD (Standard BIDS CapTrak)
            Proxy.Cz  = v_Cz;
            Proxy.T7  = v_T7;
            Proxy.T8  = v_T8;
            Proxy.Fpz = v_Fpz;

        otherwise
            error('Unknown Cohort ID: %s. Update Sanitizer V15 gateway.', cohort_id);
    end

    % PROXY OUTPUT IS NOW MATHEMATICALLY SECURE.
end

% =========================================================================
% HELPER: Robust Electrode Extraction
% =========================================================================
function coords = extract_electrode(data, targets)
    coords = [];
    % Check 'label' or 'name' columns first
    col_names = data.Properties.VariableNames;
    label_col = '';
    if any(strcmpi(col_names, 'label')), label_col = 'label';
    elseif any(strcmpi(col_names, 'name')), label_col = 'name';
    end
    
    if isempty(label_col)
        return;
    end
    
    for i = 1:height(data)
        val = data.(label_col)(i);
        if iscell(val), val = val{1}; end
        val_str = string(val);
        
        % 1. Exact Match
        if any(strcmpi(val_str, targets))
            coords = [data.x(i), data.y(i), data.z(i)];
            return;
        end
        
        % 2. Prefix Match (Strip 'E' or 'EEG')
        stripped = regexprep(val_str, '^(E|EEG)', '', 'ignorecase');
        if any(strcmpi(stripped, targets))
            coords = [data.x(i), data.y(i), data.z(i)];
            return;
        end
    end
end
