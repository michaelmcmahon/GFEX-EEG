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
% SCRIPT: Verify_12_Traps_Scoreboard.m
% Objective: Compare New Ingestion (12-Trap) vs. Old Ingestion (V13.5)
% target: Proving zero regression on working cohorts.
% =========================================================================
function Verify_12_Traps_Scoreboard()

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    raw_dir  = fullfile(root_dir, 'Data_Clean');
    proc_dir = fullfile(root_dir, 'Data_Processed'); 
    
    cohorts = {'ds007353', 'ds005811-ds005810', 'ds005795', 'ds004718', 'ds000117-ds002718', 'ds006525', 'ds007216'};
    
    fprintf('\n%s\n', repmat('=', 1, 110));
    fprintf('%-20s | %-15s | %-12s | %-12s | %-12s | %-10s\n', 'Cohort', 'Subject', 'Scale Type', 'New D3D', 'Old D3D', 'Delta (mm)');
    fprintf('%s\n', repmat('-', 1, 110));

    for c = 1:length(cohorts)
        cohort = cohorts{c};
        
        % 1. Load OLD Data
        old_file = fullfile(proc_dir, cohort, 'Master_Scalp_Proxy.mat');
        if ~exist(old_file, 'file'), fprintf('%-20s | [SKIP] No old data found.\n', cohort); continue; end
        old_data = load(old_file);
        subs = fieldnames(old_data.Master_Scalp_Proxy);
        if isempty(subs), continue; end
        sub_key = subs{1}; 
        sub_id = strrep(sub_key, '_', '-'); % Back-convert key to ID
        
        old_p = old_data.Master_Scalp_Proxy.(sub_key);
        old_cz = old_p.Cz;
        old_d3d = norm(old_p.Cz - old_p.T7) * 1000;

        % 2. Run NEW Logic (Direct from Raw)
        try
            % Label Dictionary
            dict = struct('Cz', 'Cz', 'T7', 'T7', 'T8', 'T8', 'Fpz', 'Fpz');
            if contains(cohort, 'ds000117'), dict.Cz = 'EEG034'; dict.T7 = 'EEG030'; dict.T8 = 'EEG038'; dict.Fpz = 'EEG002';
            elseif contains(cohort, 'ds006525'), dict.Cz = '80'; dict.T7 = '44'; dict.T8 = '77'; dict.Fpz = '11'; end

            % Find raw file
            deriv_root = fullfile(raw_dir, cohort, 'derivatives', 'intensity_normalization');
            if ~exist(deriv_root, 'dir'), deriv_root = fullfile(raw_dir, cohort); end
            d_t = dir(fullfile(raw_dir, cohort, sub_id, '**', '*_electrodes.tsv'));
            if isempty(d_t), d_t = dir(fullfile(deriv_root, sub_id, '**', '*_electrodes.tsv')); end
            
            valid_idx = [];
            for i = 1:length(d_t), if ~startsWith(d_t(i).name, '._'), valid_idx = [valid_idx, i]; end; end
            e_path = fullfile(d_t(valid_idx(1)).folder, d_t(valid_idx(1)).name);
            
            data = readtable(e_path, 'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
            vNames = data.Properties.VariableNames;
            ixX = find(strcmpi(vNames, 'x'), 1); ixY = find(strcmpi(vNames, 'y'), 1); ixZ = find(strcmpi(vNames, 'z'), 1);
            ixN = find(strcmpi(vNames, 'name') | strcmpi(vNames, 'label'), 1);
            raw_all = [data{:,ixX}, data{:,ixY}, data{:,ixZ}];
            lbls = string(data{:,ixN});

            iCz = find(strcmpi(lbls, dict.Cz), 1); iT7 = find(strcmpi(lbls, dict.T7), 1); iT8 = find(strcmpi(lbls, dict.T8), 1);
            iFpz = find(strcmpi(lbls, dict.Fpz), 1);
            
            % NEW Scale logic
            d_arc = norm(raw_all(iCz,:) - raw_all(iT7,:));
            scale = 1.0; type = 'Meters';
            if d_arc > 30, scale = 0.001; type = 'Millimeters';
            elseif d_arc > 5, scale = 0.01; type = 'Centimeters';
            elseif d_arc > 0.4 && d_arc < 2.0 
                if contains(cohort, 'ds005795'), scale = 0.085 / d_arc; type = 'Sphere-R1';
                else, scale = 0.1; type = 'Decimeters'; end
            elseif d_arc > 2.0 && d_arc < 5.0, scale = 0.140 / d_arc; type = 'Sphere-PI';
            elseif d_arc < 0.25, scale = 1.0; type = 'Meters'; end
            
            xyz_m = raw_all * scale;
            
            % NEW Axis logic
            if isempty(iFpz), iFpz = iCz; end
            vAnt = xyz_m(iFpz,:) - xyz_m(iCz,:);
            [~, axMax] = max(abs(vAnt));
            if axMax == 1, xyz_ras = [-xyz_m(:,2), xyz_m(:,1), xyz_m(:,3)]; % Fix Axis
            else, xyz_ras = xyz_m; end
            
            % NEW Pivot logic
            s_pivot = (xyz_ras(iT7,:) + xyz_ras(iT8,:)) / 2;
            new_xyz = xyz_ras - s_pivot;
            
            new_cz = new_xyz(iCz,:);
            new_d3d = norm(new_xyz(iCz,:) - new_xyz(iT7,:)) * 1000;
            
            % 3. Calculate REGRESSION DELTA
            delta_cz = norm(new_cz - old_cz) * 1000;
            
            fprintf('%-20s | %-15s | %-12s | %-12.1f | %-12.1f | %-10.4f\n', cohort, sub_id, type, new_d3d, old_d3d, delta_cz);

        catch ME
            fprintf('%-20s | %-15s | FAILED: %s\n', cohort, sub_id, ME.message);
        end
    end
    fprintf('%s\n', repmat('=', 1, 110));
end
