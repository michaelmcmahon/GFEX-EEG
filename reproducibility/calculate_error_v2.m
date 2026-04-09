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

function calculate_error_v1()
    % calculate_error_v2.m (V1.20 - FULLY HARDENED)
    % Objective: Resolve RAS->ALS reflection trap, apply isotropic scaling 
    % for cap thickness, and utilize explicit repmat expansion to prevent 
    % matrix dimension crashes during translation.

    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    scripts_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile('C:\MoBI_Research\fieldtrip\fieldtrip-master\external\fileexchange'));

    % Load serialized data
    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Sensor_Clouds.mat'),         'EEG_Cloud');
    load(fullfile(res_dir, 'Master_Scalp_Proxy.mat'),    'Master_Scalp_Proxy');

    mesh = load(fullfile(scripts_dir, 'ICBM152_scalp.mat'));
    V_temp = mesh.Vertices;

    subjects = fieldnames(Algorithm_Predictions);
    Method_Error = []; Baseline_Error = [];

    fprintf('\n--- PHASE 6: FINAL ERROR ANALYSIS (V1.20) ---\n');

    for i = 1:length(subjects)
        subKey = subjects{i};
        try
            pred  = Algorithm_Predictions.(subKey);
            base  = Baseline_LPA_RPA.(subKey);
            truth = Ground_Truth_MRI.(subKey);
            cloud = EEG_Cloud.(subKey); % 64x3 matrix
            p     = Master_Scalp_Proxy.(subKey);

            % 1. INITIAL FUNCTIONAL ALIGNMENT (Rough Snap)
            [~, iTop] = max(V_temp(:,3)); [~, iFrt] = max(V_temp(:,2));
            [~, iLft] = min(V_temp(:,1)); [~, iRgt] = max(V_temp(:,1));
            P_temp_anchors = [V_temp(iTop,:); V_temp(iFrt,:); V_temp(iLft,:); V_temp(iRgt,:)];
            P_subj_anchors = [p.Cz(:)'; p.NAS(:)'; p.T7(:)'; p.T8(:)'];

            % FORBID REFLECTION to prevent RAS->ALS mirror inversion
            [~, ~, t_rough] = procrustes(P_subj_anchors, P_temp_anchors, 'scaling', true, 'reflection', false);

            % REPMAT FIX: Explicitly expand the 1x3 translation vector to match the Nx3 template mesh
            translation_expanded = repmat(t_rough.c(1,:), size(V_temp, 1), 1);
            V_temp_rough = t_rough.b * V_temp * t_rough.T + translation_expanded;

            % 2. STABLE ICP ALIGNMENT
            % With reflection forbidden and dimension errors bypassed, ICP snaps without NaN
            [R_icp, T_icp] = icp(cloud', V_temp_rough');

            % Map Template Predictions (LHJ/RHJ) into ICP-aligned CapTrak Space (Vectors are 1x3, no repmat needed here)
            predL_rough = t_rough.b * pred.LHJ(:)' * t_rough.T + t_rough.c(1,:);
            predR_rough = t_rough.b * pred.RHJ(:)' * t_rough.T + t_rough.c(1,:);

            predL_captrak = (R_icp * predL_rough' + T_icp)';
            predR_captrak = (R_icp * predR_rough' + T_icp)';

            % 3. RESOLVE THE REFERENCE FRAME DISCONNECT (CapTrak RAS vs. MEG ALS)
            Base_Fid  = [base.NAS(:)'; base.LPA(:)'; base.RPA(:)'];
            Truth_Fid = [truth.NAS(:)'; truth.LPA(:)'; truth.RPA(:)'];

            % Pure rotation/translation to match coordinate frames
            [~, ~, t_frame] = procrustes(Truth_Fid, Base_Fid, 'scaling', false, 'reflection', false);

            % Transform our predictions into the Ground Truth coordinate universe
            predL_final = predL_captrak * t_frame.T + t_frame.c(1,:);
            predR_final = predR_captrak * t_frame.T + t_frame.c(1,:);

            % Transform baseline to Ground Truth universe (Apples-to-apples)
            baseL_final = base.LPA(:)' * t_frame.T + t_frame.c(1,:);
            baseR_final = base.RPA(:)' * t_frame.T + t_frame.c(1,:);

            % 4. CALCULATE ERROR (in millimeters: * 1000 since data is natively in meters)
            err_L_method = norm(predL_final - truth.LPA(:)') * 1000;
            err_R_method = norm(predR_final - truth.RPA(:)') * 1000;
            m_err = (err_L_method + err_R_method) / 2;

            err_L_base = norm(baseL_final - truth.LPA(:)') * 1000;
            err_R_base = norm(baseR_final - truth.RPA(:)') * 1000;
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
    % With the arrays fully populated, signrank will execute properly
    [p_value, ~] = signrank(Method_Error(:), Baseline_Error(:));

    fprintf('\n=== FINAL V1.20 RESULTS (N=%d) ===\n', length(Method_Error));
    fprintf('Mean Method Error: %.2f mm (SD: %.2f)\n', mean(Method_Error), std(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm (SD: %.2f)\n', mean(Baseline_Error), std(Baseline_Error));
    fprintf('Wilcoxon p-value: %.5e\n', p_value);

end