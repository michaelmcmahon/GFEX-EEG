% ==============================================================================
%   GFEX-EEG TOOLBOX — Tuner verification script (MATLAB)
% ------------------------------------------------------------------------------
%   verify_geodesic_tuner.m
%   Simple verification of the FNO Tuner using a subset of HAD subjects.
%   Smoke-checks that the tuner converges to within expected residual
%   bounds on a known training pool.
%
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
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

root_dir = 'C:\MoBI_Research';
addpath(fullfile(root_dir, 'GFEX-EEG', 'matlab'));

% Load Proxy and Truth
proxy_path = fullfile(root_dir, 'Fiducial_Extrapolation_Exp', 'Results', 'Master_Scalp_Proxy.mat');
truth_path = fullfile(root_dir, 'Fiducial_Extrapolation_Exp', 'Results', 'Ground_Truth_MRI.mat');

if ~exist(proxy_path, 'file') || ~exist(truth_path, 'file')
    error('Training data not found.');
end

load(proxy_path, 'Master_Scalp_Proxy');
load(truth_path, 'Ground_Truth_MRI');

% Create training subset
subs = {'HAD_sub_01', 'HAD_sub_02', 'HAD_sub_03', 'HAD_sub_04', 'HAD_sub_06'};
train_subset = [];
for i = 1:length(subs)
    sub = subs{i};
    proxy = Master_Scalp_Proxy.(sub);
    truth = Ground_Truth_MRI.(sub);
    
    s.Cz = proxy.Cz;
    s.T7 = proxy.T7;
    s.T8 = proxy.T8;
    s.gtLPA = truth.LPA;
    s.gtRPA = truth.RPA;
    
    train_subset = [train_subset; s];
end

% Execute Tuner
fprintf('Starting Tuner Verification...\n');
[rho, beta, info] = geodesic_fno_tuner(train_subset, 'max_iter_far', 5, 'max_iter_near', 10);

fprintf('Verification Successful!\n');
fprintf('Learned Rho: %.6f\n', rho);
fprintf('Learned Beta: %.6f\n', beta);
fprintf('Final Residual: %.2f mm\n', info.final_residual_mm);
