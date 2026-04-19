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

function [EEG, com] = pop_geodesic_rescue(EEG, rho, beta, cohort)
% POP_GEODESIC_RESCUE (V21.7 - Invincible EEGLAB Wrapper)
% Includes: Empirical Scale Shield and RAS-to-ALS Axis Transposition.
%
% Weight precedence: explicit rho/beta > cohort preset (from weight_zoo.json)
% > hard default (LEMON-tuned, 2026-04-19).

    com = '';
    if exist('geodesic_rescue', 'file') == 0
        base_path = fileparts(mfilename('fullpath'));
        addpath(fullfile(base_path, '..', '..', 'core'));
    end

    if nargin < 4, cohort = ''; end
    if nargin < 3 || isempty(beta), beta = []; end
    if nargin < 2 || isempty(rho),  rho  = []; end

    % Resolve weights: explicit > cohort preset > hard default
    if ~isempty(cohort)
        preset = load_cohort_preset(cohort);
        if isempty(rho),  rho  = preset.rho;        end
        if isempty(beta), beta = preset.beta;       end
        D_standard = preset.D_standard;
    else
        if isempty(rho),  rho  = 0.248383;  end  % Hard default (LEMON-tuned)
        if isempty(beta), beta = 0.235926; end
        D_standard = 0.1388;
    end

    % 1. Hardware-Agnostic Channel Extraction
    idx_Cz = []; idx_T7 = []; idx_T8 = [];
    for i = 1:length(EEG.chanlocs)
        lbl = EEG.chanlocs(i).labels;
        if strcmpi(lbl, 'Cz') || strcmpi(lbl, 'E36') || strcmpi(lbl, '80'), idx_Cz = i;
        elseif strcmpi(lbl, 'T7') || strcmpi(lbl, 'T3') || strcmpi(lbl, 'E45') || strcmpi(lbl, '44'), idx_T7 = i;
        elseif strcmpi(lbl, 'T8') || strcmpi(lbl, 'T4') || strcmpi(lbl, 'E108') || strcmpi(lbl, '102'), idx_T8 = i;
        end
    end
    
    if isempty(idx_Cz) || isempty(idx_T7) || isempty(idx_T8)
        error('GeodesicRescue:Nomenclature', 'Could not identify Cz, T7, or T8 anchors.');
    end

    % 2. Absolute Unit Forcing (The V21.8 Hyper-Scale Empirical Shield)
    Cz_raw = [double(EEG.chanlocs(idx_Cz).X), double(EEG.chanlocs(idx_Cz).Y), double(EEG.chanlocs(idx_Cz).Z)];
    T7_raw = [double(EEG.chanlocs(idx_T7).X), double(EEG.chanlocs(idx_T7).Y), double(EEG.chanlocs(idx_T7).Z)];
    T8_raw = [double(EEG.chanlocs(idx_T8).X), double(EEG.chanlocs(idx_T8).Y), double(EEG.chanlocs(idx_T8).Z)];

    emp_rad = norm(Cz_raw - T7_raw);
    
    if emp_rad < 1.5 
        scale_factor = 1.0;          % Natively in Meters (e.g., HAD/Boston)
    elseif emp_rad >= 1.5 && emp_rad < 30.0
        scale_factor = 0.01;         % Natively in Centimeters
    elseif emp_rad >= 30.0 && emp_rad < 500.0
        scale_factor = 0.001;        % Natively in Millimeters (e.g., Raw HK)
    else
        % THE HYPER-SCALE CATCH (Toolbox Double-Scaling Trap)
        % e.g., radius = 84,500 units. Return to ~0.0845 Meters.
        scale_factor = 0.000001; 
    end

    % Force to absolute Meter-Space for Engine
    Cz_m = Cz_raw * scale_factor;
    T7_m = T7_raw * scale_factor;
    T8_m = T8_raw * scale_factor;

    % 3. Execute Core Engine (Meter-Space Output)
    [pLHJ_m, pRHJ_m, ~] = geodesic_rescue(Cz_m, T7_m, T8_m, ...
        'rho', rho, 'beta', beta, 'D_standard', D_standard);

    % 4. Axis Transposition (RAS to ALS) and Scale Reversal
    % MNI (RAS)    = [Right, Anterior, Superior]
    % EEGLAB (ALS) = [Anterior, Left, Superior]
    % Mapping: ALS_X = MNI_Y | ALS_Y = -MNI_X | ALS_Z = MNI_Z
    LHJ_als = [pLHJ_m(2), -pLHJ_m(1), pLHJ_m(3)] / scale_factor;
    RHJ_als = [pRHJ_m(2), -pRHJ_m(1), pRHJ_m(3)] / scale_factor;

    % 5. Inject Rescued Fiducials Back into EEG.chanlocs
    num_chans = length(EEG.chanlocs);
    EEG.chanlocs(num_chans + 1).labels = 'LPA';
    EEG.chanlocs(num_chans + 1).X = LHJ_als(1);
    EEG.chanlocs(num_chans + 1).Y = LHJ_als(2);
    EEG.chanlocs(num_chans + 1).Z = LHJ_als(3);
    EEG.chanlocs(num_chans + 1).type = 'FID';

    EEG.chanlocs(num_chans + 2).labels = 'RPA';
    EEG.chanlocs(num_chans + 2).X = RHJ_als(1);
    EEG.chanlocs(num_chans + 2).Y = RHJ_als(2);
    EEG.chanlocs(num_chans + 2).Z = RHJ_als(3);
    EEG.chanlocs(num_chans + 2).type = 'FID';

    if ~isempty(cohort)
        com = sprintf('EEG = pop_geodesic_rescue(EEG, %f, %f, ''%s'');', rho, beta, cohort);
    else
        com = sprintf('EEG = pop_geodesic_rescue(EEG, %f, %f);', rho, beta);
    end
    fprintf('SUCCESS: Rescued LPA/RPA fiducials appended (Scale Factor: %.4f).\n', scale_factor);
end

