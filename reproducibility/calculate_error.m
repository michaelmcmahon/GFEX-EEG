function calculate_error()
%CALCULATE_ERROR  Per-subject HTJ error against MRI ground truth.
%
%   Loads Algorithm_Predictions, Baseline_LPA_RPA, Ground_Truth_MRI,
%   and Transformation_Matrices from the Results folder and computes
%   the per-subject HTJ error (predicted vs. MRI ground-truth) for
%   each method on the scoreboard.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Per-subject HTJ error calculator (MATLAB)
% ------------------------------------------------------------------------------
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
    root_dir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp';
    res_dir  = fullfile(root_dir, 'Results');
    
    load(fullfile(res_dir, 'Algorithm_Predictions.mat'), 'Algorithm_Predictions');
    load(fullfile(res_dir, 'Baseline_LPA_RPA.mat'),      'Baseline_LPA_RPA');
    load(fullfile(res_dir, 'Ground_Truth_MRI.mat'),      'Ground_Truth_MRI');
    load(fullfile(res_dir, 'Transformation_Matrices.mat'), 'Trans_Mat');
    
    subjects = fieldnames(Algorithm_Predictions);
    num_subj = length(subjects);
    
    Method_Error = [];
    Baseline_Error = [];
    Results_Table = {};
    
    fprintf('\n--- PHASE 6.1: FINAL SPATIAL ERROR ANALYSIS (DIMENSION FIXED) ---\n');
    
    for i = 1:num_subj
        subKey = subjects{i};
        if ~isfield(Trans_Mat, subKey), continue; end
        T = Trans_Mat.(subKey);
        
        pred = Algorithm_Predictions.(subKey);
        base = Baseline_LPA_RPA.(subKey);
        truth = Ground_Truth_MRI.(subKey);
        
        if isempty(pred.LHJ) || isempty(base.LPA) || isempty(truth.LPA), continue; end
        
        % --- THE DIMENSION-SAFE FIX: Force everything to COLUMN VECTORS (3x1) ---
        pL = pred.LHJ(:); pR = pred.RHJ(:);
        bL = base.LPA(:); bR = base.RPA(:);
        tL = truth.LPA(:); tR = truth.RPA(:);
        
        % Apply 4x4 Trans (T * [x;y;z;1])
        pL_tr = T * [pL; 1]; pR_tr = T * [pR; 1];
        bL_tr = T * [bL; 1]; bR_tr = T * [bR; 1];
        
        % Calculate L2 Norm (Distance in mm)
        % Using norm(v1 - v2) where both are 3x1 columns
        errL_method = norm(pL_tr(1:3) - tL) * 1000;
        errR_method = norm(pR_tr(1:3) - tR) * 1000;
        
        errL_base = norm(bL_tr(1:3) - tL) * 1000;
        errR_base = norm(bR_tr(1:3) - tR) * 1000;
        
        m_err = (errL_method + errR_method) / 2;
        b_err = (errL_base + errR_base) / 2;
        
        Method_Error(end+1) = m_err;
        Baseline_Error(end+1) = b_err;
        Results_Table(end+1,:) = {subKey, m_err, b_err};
    end
    
    if isempty(Method_Error)
        fprintf('   [FAIL] No subjects matched.\n'); return;
    end
    
    [p_value, h] = signrank(Method_Error(:), Baseline_Error(:));
    
    fprintf('\n=== FINAL EXPERIMENTAL RESULTS (N=%d) ===\n', length(Method_Error));
    fprintf('Mean Method Error:   %.2f mm (SD: %.2f)\n', mean(Method_Error), std(Method_Error));
    fprintf('Mean Baseline Error: %.2f mm (SD: %.2f)\n', mean(Baseline_Error), std(Baseline_Error));
    fprintf('Wilcoxon p-value:    %.5e\n', p_value);
    
    % Save data
    save(fullfile(res_dir, 'Final_Statistical_Results.mat'), 'Method_Error', 'Baseline_Error', 'p_value');
end
