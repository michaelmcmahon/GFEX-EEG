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
% calculate_error_v8.m (V1.25 - CROSS-SESSION KEY MATCHER)
% Objective: Resolves silent BIDS session key mismatches (e.g., ses-ImageNet01 
% vs. ses-MRI) by deploying regex-based fuzzy key matching across structs.
% Retains Phase 8 Anatomical Vector Normalization and Radial Telescoping.
% =========================================================================
function calculate_error_v8()

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    scripts_dir = fileparts(mfilename('fullpath'));
    
    % Add FieldTrip (or equivalent toolbox) path for ICP and other dependencies
    addpath(fullfile('C:\MoBI_Research\fieldtrip\fieldtrip-master\external\fileexchange'));

    % Load serialized data
    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Sensor_Clouds.mat'),         'EEG_Cloud');
    load(fullfile(res_dir, 'Master_Scalp_Proxy.mat'),    'Master_Scalp_Proxy');

    % Load ICBM152 Manifold Mesh
    mesh = load(fullfile(scripts_dir, 'ICBM152_scalp.mat'));
    V_temp = mesh.Vertices;

    % Extract available fieldnames for cross-referencing
    subjects     = fieldnames(Algorithm_Predictions);
    base_fields  = fieldnames(Baseline_LPA_RPA);
    truth_fields = fieldnames(Ground_Truth_MRI);
    cloud_fields = fieldnames(EEG_Cloud);
    proxy_fields = fieldnames(Master_Scalp_Proxy);

    Method_Error = []; Baseline_Error = [];

    fprintf('\n--- PHASE 8: FINAL ERROR ANALYSIS (V1.25 CROSS-SESSION MATCHING) ---\n');

    for i = 1:length(subjects)
        subKey = subjects{i};
        try
            % CRITICAL FIX: Extract the core subject ID (e.g., 'sub_01') to bypass session suffix mismatches
            baseSubID = regexp(subKey, 'sub_[\w\d]+', 'match', 'once');
            if isempty(baseSubID)
                baseSubID = subKey; % Fallback if regex misses
            end

            % Dynamically find the matching fieldnames in the other structs
            idx_base  = find(contains(base_fields, baseSubID), 1);
            idx_truth = find(contains(truth_fields, baseSubID), 1);
            idx_cloud = find(contains(cloud_fields, baseSubID), 1);
            idx_proxy = find(contains(proxy_fields, baseSubID), 1);

            % Ensure all required data exists for this subject across sessions
            if isempty(idx_base) || isempty(idx_truth) || isempty(idx_cloud) || isempty(idx_proxy)
                fprintf('Skipping %s: Missing matched data across BIDS sessions.\n', baseSubID);
                continue;
            end

            % Load the properly matched session data
            pred  = Algorithm_Predictions.(subKey);
            base  = Baseline_LPA_RPA.(base_fields{idx_base});
            truth = Ground_Truth_MRI.(truth_fields{idx_truth});
            cloud = EEG_Cloud.(cloud_fields{idx_cloud}); 
            p     = Master_Scalp_Proxy.(proxy_fields{idx_proxy});

            % 1. INITIAL FUNCTIONAL ALIGNMENT (Triad-to-Triad Rough Snap)
            P_temp_anchors = [V_temp(pred.iCz,:); V_temp(pred.iT7,:); V_temp(pred.iT8,:)];
            P_subj_anchors = [p.Cz(:)'; p.T7(:)'; p.T8(:)'];

            % Forbid reflection to prevent inside-out flipping, allow isotropic scaling
            [~, ~, t_rough] = procrustes(P_subj_anchors, P_temp_anchors, 'scaling', true, 'reflection', false);

            % Explicitly expand the 1x3 translation vector
            translation_expanded = repmat(t_rough.c(1,:), size(V_temp, 1), 1);
            V_temp_rough = t_rough.b * V_temp * t_rough.T + translation_expanded;

            % 2. STABLE ICP ALIGNMENT
            [R_icp, T_icp] = icp(cloud', V_temp_rough');

            % Map Template Predictions (LHJ/RHJ) into ICP-aligned CapTrak Space
            predL_rough = t_rough.b * pred.LHJ(:)' * t_rough.T + t_rough.c(1,:);
            predR_rough = t_rough.b * pred.RHJ(:)' * t_rough.T + t_rough.c(1,:);

            predL_captrak = (R_icp * predL_rough' + T_icp)';
            predR_captrak = (R_icp * predR_rough' + T_icp)';

            % 3. PHASE 8: ANATOMICAL VECTOR NORMALIZATION
            % 3A. Construct the Relative Vectors in CapTrak Space (Origin = NAS)
            vec_base_L = base.LPA(:)' - base.NAS(:)';
            vec_base_R = base.RPA(:)' - base.NAS(:)';

            vec_pred_L = predL_captrak - base.NAS(:)';
            vec_pred_R = predR_captrak - base.NAS(:)';

            % 3B. Construct the Relative Vectors in True MRI Space (Origin = NAS)
            vec_true_L = truth.LPA(:)' - truth.NAS(:)';
            vec_true_R = truth.RPA(:)' - truth.NAS(:)';

            % 3C. Radial Telescoping Factor
            scale_factor = norm(truth.NAS(:)') / norm(base.NAS(:)');

            vec_pred_L_scaled = vec_pred_L * scale_factor;
            vec_pred_R_scaled = vec_pred_R * scale_factor;

            vec_base_L_scaled = vec_base_L * scale_factor;
            vec_base_R_scaled = vec_base_R * scale_factor;

            % 4. CALCULATE TRUE RELATIVE SPATIAL ERROR (in millimeters: * 1000)
            err_L_method = norm(vec_pred_L_scaled - vec_true_L) * 1000;
            err_R_method = norm(vec_pred_R_scaled - vec_true_R) * 1000;
            m_err = (err_L_method + err_R_method) / 2;

            err_L_base = norm(vec_base_L_scaled - vec_true_L) * 1000;
            err_R_base = norm(vec_base_R_scaled - vec_true_R) * 1000;
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

    % 5. STATISTICAL ANALYSIS
    if isempty(Method_Error)
        fprintf('\n[CRITICAL FAILURE] Method_Error array is empty. No subjects were processed.\n');
        return;
    end
    
    [p_value, ~] = signrank(Method_Error(:), Baseline_Error(:));

    fprintf('\n=== FINAL V1.25 RESULTS (N=%d) ===\n', length(Method_Error));
    fprintf('Mean Method Error: %.2f mm (SD: %.2f)\n', mean(Method_Error), std(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm (SD: %.2f)\n', mean(Baseline_Error), std(Baseline_Error));
    fprintf('Wilcoxon p-value: %.5e\n', p_value);

end