% /*******************************************************************************
% * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
% * Version: 1.1.0-dev
% * Repository: https://github.com/michaelmcmahon/GFEX-EEG
% * License:  MIT License
% * Authors: Michael McMahon / University of Galway
% *******************************************************************************/

function [pL_corr, pR_corr, info] = apply_mlp_correction(Cz, T7, T8, pL_geo, pR_geo, mlp_path)
% APPLY_MLP_CORRECTION  Residual-learning MLP correction on top of geodesic prediction.
%
%   [pL_corr, pR_corr, info] = apply_mlp_correction(Cz, T7, T8, pL_geo, pR_geo, mlp_path)
%
%   Loads a small 1-hidden-layer MLP from `mlp_path` (a .mat file produced by
%   c_Train_LEMON_MLP.py) and returns:
%     pL_corr = pL_geo + dL      (corrected left HTJ, metres)
%     pR_corr = pR_geo + dR      (corrected right HTJ, metres)
%     info    = struct('dL', dL, 'dR', dR, 'magnitude_mm', ...)
%
%   Inputs are all 1x3 row vectors in RAS metres.
%
%   Weights file is cached between calls to avoid repeated disk IO.

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
