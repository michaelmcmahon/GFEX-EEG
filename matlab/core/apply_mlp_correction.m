function [pL_corr, pR_corr, info] = apply_mlp_correction(Cz, T7, T8, pL_geo, pR_geo, mlp_path)
%APPLY_MLP_CORRECTION  Tier 1.5 residual-correction MLP on top of geodesic prediction.
%
%   [pL_corr, pR_corr, info] = APPLY_MLP_CORRECTION(Cz, T7, T8, pL_geo, pR_geo, mlp_path)
%
%   Loads a small 1-hidden-layer MLP from MLP_PATH (a .mat file produced
%   by c_Train_LEMON_MLP.py / per-cohort training analog) and returns:
%       pL_corr = pL_geo + dL    (corrected left  HTJ, metres)
%       pR_corr = pR_geo + dR    (corrected right HTJ, metres)
%       info    = struct('dL', dL, 'dR', dR, 'magnitude_mm', ...)
%
%   Architecture: 15-input -> 64 ReLU -> 6-output (1,414 parameters).
%   Inputs are all 1x3 row vectors in RAS metres. Weights file is cached
%   between calls to avoid repeated disk IO. Inference is pure matrix
%   multiplication (no Statistics & ML Toolbox dependency).
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Tier 1.5 residual-correction MLP (MATLAB)
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

    persistent CACHE
    if isempty(CACHE) || ~isfield(CACHE, 'path') || ~strcmp(CACHE.path, mlp_path)
        if ~isfile(mlp_path)
            error('GeodesicRescue:MLPMissing', 'MLP weights file not found: %s', mlp_path);
        end
        M = load(mlp_path);
        CACHE.path   = mlp_path;
        CACHE.W1     = double(M.W1);
        CACHE.b1     = double(M.b1(:)');   % ensure 1xN row
        CACHE.W2     = double(M.W2);
        CACHE.b2     = double(M.b2(:)');
        CACHE.x_mean = double(M.x_mean(:)');
        CACHE.x_std  = double(M.x_std(:)');
        CACHE.y_mean = double(M.y_mean(:)');
        CACHE.y_std  = double(M.y_std(:)');
    end

    % 1. Build input feature vector
    x = [Cz(:)' T7(:)' T8(:)' pL_geo(:)' pR_geo(:)'];   % 1x15

    % 2. Normalise input
    x_norm = (x - CACHE.x_mean) ./ CACHE.x_std;

    % 3. Hidden layer (ReLU)
    h = x_norm * CACHE.W1 + CACHE.b1;
    h(h < 0) = 0;

    % 4. Output layer (linear)
    y_norm = h * CACHE.W2 + CACHE.b2;

    % 5. Denormalise output (residual in metres)
    y = y_norm .* CACHE.y_std + CACHE.y_mean;

    dL = y(1:3);
    dR = y(4:6);

    pL_corr = pL_geo(:)' + dL;
    pR_corr = pR_geo(:)' + dR;

    info.dL = dL;
    info.dR = dR;
    info.dL_magnitude_mm = norm(dL) * 1000;
    info.dR_magnitude_mm = norm(dR) * 1000;
end
