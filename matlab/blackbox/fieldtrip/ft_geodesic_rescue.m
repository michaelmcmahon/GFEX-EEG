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

function [elec] = ft_geodesic_rescue(cfg, elec)
% FT_GEODESIC_RESCUE (V21.8 - Invincible FieldTrip Wrapper)
% Includes: Hyper-Scale Empirical Shield and RAS-to-ALS Axis Transposition.

    if nargin < 2
        help ft_geodesic_rescue;
        return;
    end

    % Standard defaults
    rho  = ft_getopt(cfg, 'rho', 0.0054);
    beta = ft_getopt(cfg, 'beta', 1.1885);

    % 1. Hardware-Agnostic Channel Extraction
    idx_Cz = find(strcmpi(elec.label, 'Cz') | strcmpi(elec.label, 'E36') | strcmpi(elec.label, '80'));
    idx_T7 = find(strcmpi(elec.label, 'T7') | strcmpi(elec.label, 'T3') | strcmpi(elec.label, 'E45') | strcmpi(elec.label, '44'));
    idx_T8 = find(strcmpi(elec.label, 'T8') | strcmpi(elec.label, 'T4') | strcmpi(elec.label, 'E108') | strcmpi(elec.label, '102'));

    if isempty(idx_Cz) || isempty(idx_T7) || isempty(idx_T8)
        error('FT_GEODESIC_RESCUE: Nomenclature mismatch. Could not find Cz, T7, or T8.');
    end

    % 2. Absolute Unit Forcing (The V21.8 Hyper-Scale Empirical Shield)
    Cz_raw = elec.elecpos(idx_Cz(1), :);
    T7_raw = elec.elecpos(idx_T7(1), :);
    T8_raw = elec.elecpos(idx_T8(1), :);

    emp_rad = norm(Cz_raw - T7_raw);
    if emp_rad < 1.5 
        scale_f = 1.0;
    elseif emp_rad >= 1.5 && emp_rad < 30.0
        scale_f = 0.01;
    elseif emp_rad >= 30.0 && emp_rad < 500.0
        scale_f = 0.001;
    else
        % THE HYPER-SCALE CATCH (Toolbox Double-Scaling)
        scale_f = 0.000001;
    end

    % Force to absolute Meter-Space
    Cz_m = Cz_raw * scale_f; T7_m = T7_raw * scale_f; T8_m = T8_raw * scale_f;

    % 3. Execute Geodesic Core Engine 
    [pLHJ_m, pRHJ_m, ~] = geodesic_rescue(Cz_m, T7_m, T8_m, 'rho', rho, 'beta', beta);

    % 4. Axis Transposition (RAS to ALS) and Scale Reversal
    % MNI (RAS)    = [Right, Anterior, Superior]
    % FieldTrip (ALS) = [Anterior, Left, Superior]
    LHJ_als = [pLHJ_m(2), -pLHJ_m(1), pLHJ_m(3)] / scale_f;
    RHJ_als = [pRHJ_m(2), -pRHJ_m(1), pRHJ_m(3)] / scale_f;

    % 5. Inject into FieldTrip arrays
    elec.elecpos(end+1, :) = LHJ_als;
    elec.chanpos(end+1, :) = LHJ_als;
    elec.label{end+1}      = 'LPA';

    elec.elecpos(end+1, :) = RHJ_als;
    elec.chanpos(end+1, :) = RHJ_als;
    elec.label{end+1}      = 'RPA';
    
    fprintf('FT_GEODESIC_RESCUE: SUCCESS (Empirical Scale Factor: %.6f).\n', scale_f);
end

