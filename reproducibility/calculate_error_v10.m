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
% calculate_error_v10.m (V1.27 - INTRINSIC ANATOMICAL DEFLATION)
% Objective: Resolves the physical Cap Thickness Standoff (~18mm). Utilizes 
% the inverse of the Procrustes template scalar (1/t_rough.b) to deflate 
% CapTrak predictions back to bare-scalp anatomy without MRI data leakage.
% =========================================================================
function calculate_error_v10()

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    scripts_dir = fileparts(mfilename('fullpath'));
    
    addpath(fullfile('C:\MoBI_Research\fieldtrip\fieldtrip-master\external\fileexchange'));

    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Sensor_Clouds.mat'),         'EEG_Cloud');
    load(fullfile(res_dir, 'Master_Scalp_Proxy.mat'),    'Master_Scalp_Proxy');

    mesh = load(fullfile(scripts_dir, 'ICBM152_scalp.mat'));
    V_temp = mesh.Vertices;

    subjects     = fieldnames(Algorithm_Predictions);
    base_fields  = fieldnames(Baseline_LPA_RPA);
    truth_fields = fieldnames(Ground_Truth_MRI);
    cloud_fields = fieldnames(EEG_Cloud);
    proxy_fields = fieldnames(Master_Scalp_Proxy);

    Method_Error = []; Baseline_Error = [];

    fprintf('\n--- PHASE 9: FINAL ERROR ANALYSIS (V1.27 INTRINSIC DEFLATION) ---\n');

    for i = 1:length(subjects)
        subKey = subjects{i};
        try
            baseSubID = regexp(subKey, 'sub_[\w\d]+', 'match', 'once');
            if isempty(baseSubID), baseSubID = subKey; end

            idx_base  = find(contains(base_fields, baseSubID), 1);
            idx_truth = find(contains(truth_fields, baseSubID), 1);
            idx_cloud = find(contains(cloud_fields, baseSubID), 1);
            idx_proxy = find(contains(proxy_fields, baseSubID), 1);

            if isempty(idx_base) || isempty(idx_truth) || isempty(idx_cloud) || isempty(idx_proxy)
                continue;
            end

            pred  = Algorithm_Predictions.(subKey);
            base  = Baseline_LPA_RPA.(base_fields{idx_base});
            truth = Ground_Truth_MRI.(truth_fields{idx_truth});
            cloud = EEG_Cloud.(cloud_fields{idx_cloud}); 
            p     = Master_Scalp_Proxy.(proxy_fields{idx_proxy});

            % 1. INITIAL FUNCTIONAL ALIGNMENT
            P_temp_anchors = [V_temp(pred.iCz,:); V_temp(pred.iT7,:); V_temp(pred.iT8,:)];
            P_subj_anchors = [p.Cz(:)'; p.T7(:)'; p.T8(:)'];

            [~, ~, t_rough] = procrustes(P_subj_anchors, P_temp_anchors, 'scaling', true, 'reflection', false);
            translation_expanded = repmat(t_rough.c(1,:), size(V_temp, 1), 1);
            V_temp_rough = t_rough.b * V_temp * t_rough.T + translation_expanded;

            % 2. STABLE ICP ALIGNMENT
            [R_icp, T_icp] = icp(cloud', V_temp_rough');

            predL_rough = t_rough.b * pred.LHJ(:)' * t_rough.T + t_rough.c(1,:);
            predR_rough = t_rough.b * pred.RHJ(:)' * t_rough.T + t_rough.c(1,:);

            predL_captrak = (R_icp * predL_rough' + T_icp)';
            predR_captrak = (R_icp * predR_rough' + T_icp)';

            % 3. PHASE 8: BILATERAL ORIGIN NORMALIZATION
            pred_Origin = (predL_captrak + predR_captrak) / 2;
            true_Origin = (truth.LPA(:)' + truth.RPA(:)') / 2;
            base_Origin = (base.LPA(:)'  + base.RPA(:)') / 2;

            predL_rel = predL_captrak - pred_Origin;
            predR_rel = predR_captrak - pred_Origin;

            trueL_rel = truth.LPA(:)' - true_Origin;
            trueR_rel = truth.RPA(:)' - true_Origin;

            baseL_rel = base.LPA(:)' - base_Origin;
            baseR_rel = base.RPA(:)' - base_Origin;

            % 4. PHASE 9: INTRINSIC ANATOMICAL DEFLATION
            % t_rough.b is the exact scalar that inflated the bare-scalp ICBM152 
            % to the CapTrak sensor cloud (cap + hair). We use its inverse to 
            % mathematically deflate the predictions back to flesh proportions.
            intrinsic_deflation = 1 / t_rough.b;

            predL_rel_deflated = predL_rel * intrinsic_deflation;
            predR_rel_deflated = predR_rel * intrinsic_deflation;
            
            % Deflate the Baseline purely to evaluate the scaled Dummy Sphere geometry
            baseL_rel_deflated = baseL_rel * intrinsic_deflation;
            baseR_rel_deflated = baseR_rel * intrinsic_deflation;

            % 5. CALCULATE TRUE RELATIVE SPATIAL ERROR (in millimeters: * 1000)
            err_L_method = norm(predL_rel_deflated - trueL_rel) * 1000;
            err_R_method = norm(predR_rel_deflated - trueR_rel) * 1000;
            m_err = (err_L_method + err_R_method) / 2;

            err_L_base = norm(baseL_rel_deflated - trueL_rel) * 1000;
            err_R_base = norm(baseR_rel_deflated - trueR_rel) * 1000;
            b_err = (err_L_base + err_R_base) / 2;

            if ~isnan(m_err) && ~isnan(b_err)
                Method_Error(end+1)   = m_err;
                Baseline_Error(end+1) = b_err;
            end
        catch ME
            fprintf('Failed subject %s: %s\n', subKey, ME.message);
            continue;
        end
    end

    if isempty(Method_Error)
        fprintf('\n[CRITICAL FAILURE] Method_Error array is empty. No subjects were processed.\n');
        return;
    end
    
    [p_value, ~] = signrank(Method_Error(:), Baseline_Error(:));

    fprintf('\n=== FINAL V1.27 RESULTS (N=%d) ===\n', length(Method_Error));
    fprintf('Mean Method Error: %.2f mm (SD: %.2f)\n', mean(Method_Error), std(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm (SD: %.2f)\n', mean(Baseline_Error), std(Baseline_Error));
    fprintf('Wilcoxon p-value: %.5e\n', p_value);

end