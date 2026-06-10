function [elec] = ft_geodesic_rescue(cfg, elec)
%FT_GEODESIC_RESCUE  FieldTrip wrapper for GFEX-EEG geodesic rescue.
%
%   elec = FT_GEODESIC_RESCUE(cfg, elec)
%
%   FieldTrip-flavoured wrapper around the GFEX-EEG core engine.
%   Generates missing LHJ / RHJ fiducial coordinates on the input
%   electrode definition `elec`, with Hyper-Scale Empirical Shield
%   and automatic RAS-to-ALS axis transposition.
%
%   Configuration fields:
%     cfg.cohort     - Cohort preset tag (see weight_zoo.json).
%     cfg.rho        - Explicit override (highest precedence).
%     cfg.beta       - Explicit override.
%     cfg.D_standard - Explicit override.
%
%   Weight precedence: explicit cfg.rho/beta/D_standard > cfg.cohort
%   preset > hard default (LEMON-tuned, 2026-04-19).
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — FieldTrip black-box wrapper (MATLAB, V21.8)
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
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
% ==============================================================================

    if nargin < 2
        help ft_geodesic_rescue;
        return;
    end

    % Weight resolution: explicit cfg.rho/beta/D_standard > cfg.cohort preset > hard default
    cohort = ft_getopt(cfg, 'cohort', '');
    rho_cfg        = ft_getopt(cfg, 'rho',        []);
    beta_cfg       = ft_getopt(cfg, 'beta',       []);
    Dstd_cfg       = ft_getopt(cfg, 'D_standard', []);

    rho        = 0.248383; % Hard default (LEMON-tuned, 2026-04-19)
    beta       = 0.235926;
    D_standard = 0.1388;
    if ~isempty(cohort)
        preset = load_cohort_preset(cohort);
        rho = preset.rho; beta = preset.beta; D_standard = preset.D_standard;
    end
    if ~isempty(rho_cfg),  rho        = rho_cfg;  end
    if ~isempty(beta_cfg), beta       = beta_cfg; end
    if ~isempty(Dstd_cfg), D_standard = Dstd_cfg; end

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
    [pLHJ_m, pRHJ_m, ~] = geodesic_rescue(Cz_m, T7_m, T8_m, ...
        'rho', rho, 'beta', beta, 'D_standard', D_standard);

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

