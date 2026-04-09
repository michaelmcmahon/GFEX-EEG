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
% calculate_error_v7.m (V1.25 - FLESH NASION VECTOR NORMALIZATION)
% Objective: Resolve the "Nasion Anchor Problem" (Cap Ring vs. Flesh Nose).
% Maps CapTrak predictions directly into True MRI Space via Mid-Ear 
% Centering and Radial Telescoping, allowing both vectors to be anchored 
% perfectly to the True Flesh Nasion.
% =========================================================================
function calculate_error_v7()

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    scripts_dir = fileparts(mfilename('fullpath'));
    
    % Add FieldTrip (or equivalent toolbox) path for ICP
    addpath(fullfile('C:\MoBI_Research\fieldtrip\fieldtrip-master\external\fileexchange'));

    % Load serialized data
    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Sensor_Clouds.mat'),         'EEG_Cloud');
    load(fullfile(res_dir, 'Master_Scalp_Proxy.mat'),    'Master_Scalp_Proxy');

    % Load ICBM152 Manifold Mesh (Used as the MRI-free Scalp Proxy target)
    mesh = load(fullfile(scripts_dir, 'ICBM152_scalp.mat'));
    V_temp = mesh.Vertices;

    subjects = fieldnames(Algorithm_Predictions);
    Method_Error = []; Baseline_Error = [];

    fprintf('\n--- PHASE 9: FINAL ERROR ANALYSIS (V1.25 FLESH NASION ANCHOR) ---\n');

    for i = 1:length(subjects)
        subKey = subjects{i};
        try
            pred  = Algorithm_Predictions.(subKey);
            base  = Baseline_LPA_RPA.(subKey);
            truth = Ground_Truth_MRI.(subKey);
            cloud = EEG_Cloud.(subKey); % 64x3 matrix
            p     = Master_Scalp_Proxy.(subKey);

            % 1. INITIAL FUNCTIONAL ALIGNMENT (Triad-to-Triad Rough Snap)
            P_temp_anchors = [V_temp(pred.iCz,:); V_temp(pred.iT7,:); V_temp(pred.iT8,:)];
            P_subj_anchors = [p.Cz(:)'; p.T7(:)'; p.T8(:)'];

            % Forbid reflection to prevent inside-out flipping, allow scaling
            [~, ~, t_rough] = procrustes(P_subj_anchors, P_temp_anchors, 'scaling', true, 'reflection', false);

            % Explicitly expand the 1x3 translation vector
            translation_expanded = repmat(t_rough.c(1,:), size(V_temp, 1), 1);
            V_temp_rough = t_rough.b * V_temp * t_rough.T + translation_expanded;

            % 2. STABLE ICP ALIGNMENT (Snap Template to EEG Cloud)
            [R_icp, T_icp] = icp(cloud', V_temp_rough');

            % Map Template Predictions (LHJ/RHJ) into ICP-aligned CapTrak Space
            predL_rough = t_rough.b * pred.LHJ(:)' * t_rough.T + t_rough.c(1,:);
            predR_rough = t_rough.b * pred.RHJ(:)' * t_rough.T + t_rough.c(1,:);

            predL_captrak = (R_icp * predL_rough' + T_icp)';
            predR_captrak = (R_icp * predR_rough' + T_icp)';

            % =====================================================================
            % 3. PHASE 9: FLESH NASION VECTOR NORMALIZATION (V1.25)
            % Objective: Map predictions into True MRI space to anchor vectors 
            % to the Flesh Nasion, bypassing the 42mm Dummy Nasion diagonal drop.
            % =====================================================================

            % 3A. Center CapTrak Predictions on their own mid-ear origin (Z=0)
            mid_pred = (predL_captrak + predR_captrak) / 2;
            predL_centered = predL_captrak - mid_pred;
            predR_centered = predR_captrak - mid_pred;
            
            % Center Baseline on its own origin
            mid_base = (base.LPA(:)' + base.RPA(:)') / 2;
            baseL_centered = base.LPA(:)' - mid_base;
            baseR_centered = base.RPA(:)' - mid_base;

            % 3B. Apply Radial Telescoping Factor to Map to MRI Scale
            scale_factor = norm(truth.NAS(:)') / norm(base.NAS(:)');
            predL_mapped = predL_centered * scale_factor;
            predR_mapped = predR_centered * scale_factor;
            
            baseL_mapped = baseL_centered * scale_factor;
            baseR_mapped = baseR_centered * scale_factor;

            % 3C. Anchor Vectors to the Aligned Flesh Nose
            % Because the coordinates are now mid-ear centered and deflated, 
            % they physically reside in the True MRI Space.
            Aligned_Flesh_Nose = truth.NAS(:)';

            vec_pred_L = predL_mapped - Aligned_Flesh_Nose;
            vec_pred_R = predR_mapped - Aligned_Flesh_Nose;

            vec_base_L = baseL_mapped - Aligned_Flesh_Nose;
            vec_base_R = baseR_mapped - Aligned_Flesh_Nose;

            % 3D. Construct True MRI Vectors (Also anchored to Flesh Nose)
            vec_true_L = truth.LPA(:)' - truth.NAS(:)';
            vec_true_R = truth.RPA(:)' - truth.NAS(:)';

            % 4. CALCULATE TRUE RELATIVE SPATIAL ERROR (in millimeters: * 1000)
            err_L_method = norm(vec_pred_L - vec_true_L) * 1000;
            err_R_method = norm(vec_pred_R - vec_true_R) * 1000;
            m_err = (err_L_method + err_R_method) / 2;

            err_L_base = norm(vec_base_L - vec_true_L) * 1000;
            err_R_base = norm(vec_base_R - vec_true_R) * 1000;
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