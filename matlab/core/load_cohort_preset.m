% /*******************************************************************************
% * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
% * Version: 1.0.0
% * Repository: https://github.com/michaelmcmahon/GFEX-EEG
% * License:  MIT License
% * Authors: Michael McMahon / University of Galway
% *******************************************************************************/

function preset = load_cohort_preset(tag, zoo_path)
% LOAD_COHORT_PRESET  Retrieve a cohort-specific FNO weight preset from the zoo.
%
%   preset = load_cohort_preset(tag)
%   preset = load_cohort_preset(tag, zoo_path)
%
%   Inputs:
%     tag       - Cohort tag string (e.g. 'LEMON_Polhemus_Adult', 'lemon',
%                 'default'). Aliases are resolved transparently.
%     zoo_path  - (Optional) Path to weight_zoo.json. Defaults to
%                 <matlab/core>/../../data/weight_zoo.json.
%
%   Output:
%     preset    - Struct with fields: rho, beta, D_standard, tag, status,
%                 description (when available).
%
%   Errors:
%     - If tag is not found, lists available production tags.
%     - If tag has status == 'pending_tune', raises with note explaining why.

    if nargin < 2 || isempty(zoo_path)
        base_path = fileparts(mfilename('fullpath'));
        zoo_path  = fullfile(base_path, '..', '..', 'data', 'weight_zoo.json');
    end

    if ~isfile(zoo_path)
        error('GeodesicRescue:ZooMissing', ...
            'weight_zoo.json not found at %s', zoo_path);
    end

    zoo = jsondecode(fileread(zoo_path));

    % Resolve aliases (case-insensitive)
    resolved = tag;
    if isfield(zoo, 'aliases')
        alias_names = fieldnames(zoo.aliases);
        for i = 1:numel(alias_names)
            if strcmpi(alias_names{i}, tag)
                resolved = zoo.aliases.(alias_names{i});
                break;
            end
        end
    end

    if ~isfield(zoo.weight_sets, resolved)
        available = fieldnames(zoo.weight_sets);
        production_tags = {};
        for i = 1:numel(available)
            entry = zoo.weight_sets.(available{i});
            if isfield(entry, 'status') && strcmp(entry.status, 'production')
                production_tags{end+1} = available{i}; %#ok<AGROW>
            end
        end
        error('GeodesicRescue:UnknownCohort', ...
            ['Cohort tag "%s" not found in weight zoo.\n' ...
             'Available production tags: %s\n' ...
             'Available aliases: see weight_zoo.json'], ...
            tag, strjoin(production_tags, ', '));
    end

    entry = zoo.weight_sets.(resolved);

    if isfield(entry, 'status') && strcmp(entry.status, 'pending_tune')
        note = '';
        if isfield(entry, 'note'), note = entry.note; end
        error('GeodesicRescue:PendingTune', ...
            ['Cohort "%s" is marked pending_tune (no production weights yet).\n' ...
             'Note: %s\n' ...
             'Use Tier 2 FNO tuner on pilot data with MRI HTJ ground truth, ' ...
             'or fall back to the default LEMON preset.'], ...
            resolved, note);
    end

    preset.tag         = resolved;
    preset.rho         = entry.rho;
    preset.beta        = entry.beta;
    preset.D_standard  = entry.D_standard;
    preset.status      = entry.status;
    if isfield(entry, 'description'), preset.description = entry.description; end
    if isfield(entry, 'mlp_weights_file'), preset.mlp_weights_file = entry.mlp_weights_file; end
    if isfield(entry, 'mlp_holdout_mean_err_mm'), preset.mlp_holdout_mean_err_mm = entry.mlp_holdout_mean_err_mm; end
end
