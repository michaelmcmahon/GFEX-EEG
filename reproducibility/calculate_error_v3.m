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
% PHASE 5: STATISTICAL VALIDATION AND SCOREBOARD GENERATOR V3.0
% Objective: Quantify optimized accuracy and verify vs. baseline.
% Logic: Euclidean L2 Norm in mm. Wilcoxon Rank-Sum for p-values.
% =========================================================================
function calculate_error_v3()

    res_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp\Results';
    scripts_dir = fileparts(mfilename('fullpath'));
    load(fullfile(res_dir, 'Master_Cohort_Training_Data.mat'), 'Master_Scalp_Proxy', 'Ground_Truth_MRI');
    load(fullfile(res_dir, 'FNO_Optimization_Final.mat'), 'optimal_theta');
    
    rho = optimal_theta(1);
    beta = optimal_theta(2);
    mesh_path = fullfile(scripts_dir, 'ICBM152_scalp.mat');
    
    % Template Constants (Meters)
    T_Cz = [0.011240, 0.025921, 0.141134];
    T_T7 = [-0.089174, -0.001327, -0.006348];
    T_T8 = [0.096880, -0.014286, -0.005819];

    subs = fieldnames(Master_Scalp_Proxy);
    Method_Errors = zeros(length(subs), 1);
    Baseline_Errors = zeros(length(subs), 1);

    fprintf('\nPhase 5: Calculating Final Error Scoreboard (N=%d)...\n', length(subs));

    for i = 1:length(subs)
        sub = subs{i};
        t = Ground_Truth_MRI.(sub);
        
        % 1. OPTIMIZED PREDICTION (FNO learned parameters)
        [Lp, Rp] = predict_helix_tragus_junctions_fno(T_Cz, T_T7, T_T8, mesh_path, rho, beta);
        
        % 2. BASELINE PREDICTION (Axis Flatline - T7/T8 as HTJ)
        % This represents the current clinical error if proxies are used as HTJs.
        L_base = T_T7; R_base = T_T8;
        
        % 3. CALCULATE EUCLIDEAN DISTANCE (mm)
        err_method = (norm(Lp - t.LPA_template) + norm(Rp - t.RPA_template)) / 2 * 1000;
        err_base   = (norm(L_base - t.LPA_template) + norm(R_base - t.RPA_template)) / 2 * 1000;
        
        Method_Errors(i) = err_method;
        Baseline_Errors(i) = err_base;
    end

    % 4. STATISTICAL VALIDATION
    p_val = signrank(Method_Errors, Baseline_Errors);
    mean_err = mean(Method_Errors);
    sd_err = std(Method_Errors);
    improvement = mean(Baseline_Errors) - mean_err;

    fprintf('\n======================================================\n');
    fprintf('  FINAL RESULTS SCOREBOARD\n');
    fprintf('  Mean Method Error:   %.2f mm (SD: %.2f)\n', mean_err, sd_err);
    fprintf('  Mean Baseline Error: %.2f mm\n', mean(Baseline_Errors));
    fprintf('  Net Improvement:     %.2f mm reduction\n', improvement);
    fprintf('  Wilcoxon p-value:    %.5e\n', p_val);
    fprintf('======================================================\n');
    
    save(fullfile(res_dir, 'Final_Statistical_Scoreboard.mat'), 'Method_Errors', 'Baseline_Errors', 'p_val');
end
