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

﻿function calculate_error_v1_DEBUG()
    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    addpath('C:\MoBI_Research\fieldtrip\fieldtrip-master\external\fileexchange');
    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Sensor_Clouds.mat'),         'EEG_Cloud');
    mesh = load('ICBM152_scalp.mat');
    subjects = fieldnames(Algorithm_Predictions);
    Method_Error = []; Baseline_Error = [];
    fprintf('\n--- PHASE 6: FINAL ERROR ANALYSIS (ICP DEBUG) ---\n');
    for i = 1:length(subjects)
        subKey = subjects{i};
        fprintf('Processing %s...\n', subKey);
        pred = Algorithm_Predictions.(subKey);
        base = Baseline_LPA_RPA.(subKey);
        truth = Ground_Truth_MRI.(subKey);
        cloud = EEG_Cloud.(subKey);
        [R, T] = icp(mesh.Vertices', cloud');
        predL = R * pred.LHJ(:) + T;
        predR = R * pred.RHJ(:) + T;
        baseL = R * base.LPA(:) + T;
        baseR = R * base.RPA(:) + T;
        m_err = (norm(predL - truth.LPA(:)) + norm(predR - truth.RPA(:))) / 2;
        b_err = (norm(baseL - truth.LPA(:)) + norm(baseR - truth.RPA(:))) / 2;
        Method_Error(end+1) = m_err;
        Baseline_Error(end+1) = b_err;
    end
    [p_value, h] = signrank(Method_Error, Baseline_Error);
    fprintf('\n=== FINAL ICP RESULTS (N=%d) ===\n', length(Method_Error));
    fprintf('Mean Method Error:   %.2f mm\n', mean(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm\n', mean(Baseline_Error));
    fprintf('Wilcoxon p-value:    %.5e\n', p_value);
end
