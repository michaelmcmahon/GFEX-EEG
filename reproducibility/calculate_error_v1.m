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

﻿function calculate_error_v1()
    % calculate_error_v1.m (V1.19 - REFLECTION HARDENED)
    % Objective: Resolve RAS->ALS axis conflict by forbidding reflection.
    
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
    
    subjects = fieldnames(Algorithm_Predictions);
    Method_Error = []; Baseline_Error = []; Results_Table = {};
    
    fprintf('\n--- PHASE 6: FINAL ERROR ANALYSIS (V1.19 REFLECTION HARDENED) ---\n');
    
    for i = 1:numel(subjects)
        subKey = subjects{i};
        try
            pred = Algorithm_Predictions.(subKey); base = Baseline_LPA_RPA.(subKey); 
            truth = Ground_Truth_MRI.(subKey); cloud = EEG_Cloud.(subKey); p = Master_Scalp_Proxy.(subKey);
            
            % 1. ROUGH ALIGNMENT (Scaling=True, Reflection=FALSE)
            % Forbids mirror-flips during RAS to ALS conversion
            [~, iTop] = max(V_temp(:,3)); [~, iFrt] = max(V_temp(:,2)); 
            [~, iLft] = min(V_temp(:,1)); [~, iRgt] = max(V_temp(:,1));
            MRI_Func = [V_temp(iTop,:); V_temp(iFrt,:); V_temp(iLft,:); V_temp(iRgt,:)];
            EEG_Func = [p.Cz(:)'; p.NAS(:)'; p.T7(:)'; p.T8(:)'];
            
            [~, ~, Trans_Pre] = procrustes(MRI_Func, EEG_Func, 'scaling', true, 'reflection', false);
            
            % Scaled and Rotated Rough Cloud
            cloud_r = Trans_Pre.b * cloud * Trans_Pre.T + repmat(Trans_Pre.c(1,:), size(cloud, 1), 1);
            
            % 2. STABLE ICP REFINEMENT
            [R_icp, T_icp] = icp(V_temp', cloud_r');
            
            % 3. TEMPLATE-TO-SUBJECT MAPPING
            P_phys_r = Trans_Pre.b * EEG_Func([1,3,4], :) * Trans_Pre.T + repmat(Trans_Pre.c(1,:), 3, 1);
            P_ali = (R_icp * P_phys_r' + T_icp)';
            P_t = [V_temp(pred.iCz,:); V_temp(pred.iT7,:); V_temp(pred.iT8,:)];
            [~, ~, t_map] = procrustes(P_ali, P_t, 'scaling', true, 'reflection', false);
            
            fL = (t_map.b * pred.LHJ(:)' * t_map.T + t_map.c(1,:))'; 
            fR = (t_map.b * pred.RHJ(:)' * t_map.T + t_map.c(1,:))';
            
            % Baseline in Aligned Space
            bL_r = (Trans_Pre.b * base.LPA(:)' * Trans_Pre.T + Trans_Pre.c(1,:))';
            bR_r = (Trans_Pre.b * base.RPA(:)' * Trans_Pre.T + Trans_Pre.c(1,:))';
            bL_ali = R_icp * bL_r + T_icp;
            bR_ali = R_icp * bR_r + T_icp;
            
            % 4. FINAL ERROR CALC (mm)
            m_err = (norm(fL - truth.LPA(:)) + norm(fR - truth.RPA(:))) / 2 * 1000;
            b_err = (norm(bL_ali - truth.LPA(:)) + norm(bR_ali - truth.RPA(:))) / 2 * 1000;
            
            if ~isnan(m_err) && ~isnan(b_err)
                Method_Error(end+1) = m_err;
                Baseline_Error(end+1) = b_err;
                Results_Table(end+1,:) = {subKey, m_err, b_err};
            end
        catch ME
            fprintf('   [FAIL] %s: %s\n', subKey, ME.message);
        end
    end
    [p_val, ~] = signrank(Method_Error(:), Baseline_Error(:));
    fprintf('\n=== FINAL V1.19 RESULTS (N=%d) ===\n', numel(Method_Error));
    fprintf('Mean Method Error:   %.2f mm (SD: %.2f)\n', mean(Method_Error), std(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm (SD: %.2f)\n', mean(Baseline_Error), std(Baseline_Error));
    fprintf('Wilcoxon p-value:    %.5e\n', p_val);
end
