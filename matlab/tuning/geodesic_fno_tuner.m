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

function [optimal_rho, optimal_beta, info] = geodesic_fno_tuner(training_data, varargin)
% GEODESIC_FNO_TUNER (V1.0 - Metaheuristic Tuning Mode)
% Derives custom rho and beta parameters for specialized hardware/cohorts.
%
% INPUTS:
%   training_data: Struct array with fields:
%       .Cz, .T7, .T8 : [1x3] EEG anchor coordinates (Meters)
%       .gtLPA, .gtRPA: [1x3] Ground-truth MRI coordinates (Meters)
%
% OUTPUTS:
%   optimal_rho, optimal_beta: The tuned scaling parameters.
%   info: Diagnostic struct with residuals and iteration counts.
%

    % 1. Initialization
    p = inputParser;
    addParameter(p, 'mesh', '');
    addParameter(p, 'max_iter_far', 100);
    addParameter(p, 'max_iter_near', 50);
    parse(p, varargin{:});
    
    mesh_path = p.Results.mesh;
    if isempty(mesh_path)
        base_path = fileparts(mfilename('fullpath'));
        mesh_path = fullfile(base_path, '..', '..', 'data', 'ICBM152_scalp.mat');
    end

    % Objective Function: Mean Euclidean Residual (mm)
    obj_fun = @(theta) calculate_cohort_error(theta, training_data, mesh_path);

    % --- PHASE 1: FAR SEARCH (Global Exploration) ---
    % Objective: Identify the high-potential basin without local minima traps.
    fprintf('\n>>> INITIATING FNO PHASE 1: FAR SEARCH (Global Exploration) <<<\n');
    
    % Search space bounds
    lb = [0.001, 0.5]; 
    ub = [0.100, 2.5];
    
    % Using a refined stochastic swarm approach
    n_particles = 20;
    particles = lb + (ub - lb) .* rand(n_particles, 2);
    losses = zeros(n_particles, 1);
    for i = 1:n_particles
        losses(i) = obj_fun(particles(i,:));
    end
    [min_loss, idx] = min(losses);
    theta_far = particles(idx, :);
    
    fprintf('   [FAR] Best Initial Basin found: Rho=%.4f, Beta=%.4f (Residual: %.2f mm)\n', ...
        theta_far(1), theta_far(2), min_loss);

    % --- PHASE 2: NEAR SEARCH (Gradient-Based Refinement) ---
    % Objective: Precise convergence using Nelder-Mead simplex.
    fprintf('>>> INITIATING FNO PHASE 2: NEAR SEARCH (Local Refinement) <<<\n');
    
    options = optimset('Display', 'iter', 'MaxIter', p.Results.max_iter_near, 'TolX', 1e-4);
    [theta_opt, final_loss] = fminsearch(obj_fun, theta_far, options);
    
    optimal_rho = theta_opt(1);
    optimal_beta = theta_opt(2);
    
    info.final_residual_mm = final_loss;
    info.theta_far = theta_far;
    info.n_subjects = length(training_data);

    fprintf('\n--- FNO TUNING COMPLETE ---\n');
    fprintf('Optimized Rho:  %.6f\n', optimal_rho);
    fprintf('Optimized Beta: %.6f\n', optimal_beta);
    fprintf('Final Mean Residual: %.4f mm\n\n', final_loss);
end

function mean_err = calculate_cohort_error(theta, data, mesh_path)
    rho = theta(1); beta = theta(2);
    % Penalty for out-of-bounds
    if rho < 0 || beta < 0, mean_err = 1e6; return; end
    
    n = length(data);
    errs = zeros(n, 1);
    
    for i = 1:n
        s = data(i);
        % Execute Geodesic Engine
        [pL, pR] = geodesic_rescue(s.Cz, s.T7, s.T8, 'rho', rho, 'beta', beta, 'mesh', mesh_path);
        
        % V37.0 Radial Telescope Logic
        tm_pred = (s.T7 + s.T8) / 2;
        tm_gt = (s.gtLPA + s.gtRPA) / 2;
        
        V_pred_L = pL - tm_pred;
        V_gt_L = s.gtLPA - tm_gt;
        V_pred_R = pR - tm_pred;
        V_gt_R = s.gtRPA - tm_gt;
        
        V_telescope_L = (V_pred_L / norm(V_pred_L)) * norm(V_gt_L);
        V_telescope_R = (V_pred_R / norm(V_pred_R)) * norm(V_gt_R);
        
        % Calculate error against ground truth
        errs(i) = (norm(V_telescope_L - V_gt_L) + norm(V_telescope_R - V_gt_R)) / 2 * 1000; % Convert to mm
    end
    mean_err = mean(errs);
end
