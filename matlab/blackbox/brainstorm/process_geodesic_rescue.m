function varargout = process_geodesic_rescue( varargin )
%PROCESS_GEODESIC_RESCUE  Brainstorm process node wrapper for GFEX-EEG geodesic rescue.
%
%   Brainstorm "Custom > Anatomy" process node that runs the GFEX-EEG
%   geodesic walk on the active channel file and writes the predicted
%   LHJ / RHJ coordinates back into Brainstorm's SCS frame.
%
%   Options exposed in the Brainstorm GUI: cohort preset, explicit
%   rho / beta override, optional Tier 1.5 MLP residual correction.
%   Weight precedence (highest first): explicit rho/beta > cohort
%   preset (from weight_zoo.json) > hard default (LEMON-tuned).
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Brainstorm black-box wrapper (MATLAB, V21.8)
% ------------------------------------------------------------------------------
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
%   Issues     : https://github.com/michaelmcmahon/GFEX-EEG/issues
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

    eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Geodesic Origin Rescue (V21.8)';
    sProcess.Category    = 'Custom';
    sProcess.SubCategory = 'Anatomy';
    sProcess.Index       = 1000;
    sProcess.InputTypes  = {'channel'};
    sProcess.OutputTypes = {'channel'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.cohort.Comment = 'Cohort preset (empty = use Rho/Beta below):';
    sProcess.options.cohort.Type    = 'text';
    sProcess.options.cohort.Value   = '';
    sProcess.options.rho.Comment = 'Optimal Ratio (Rho):';
    sProcess.options.rho.Type    = 'value';
    sProcess.options.rho.Value   = 0.248383;
    sProcess.options.beta.Comment = 'Scaling Limit (Beta):';
    sProcess.options.beta.Type    = 'value';
    sProcess.options.beta.Value   = 0.235926;
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {sInputs.FileName};
    ChannelMat = in_bst_channel(sInputs.FileName);
    
    % 1. Hardware-Agnostic Channel Extraction
    iCz = find(strcmpi({ChannelMat.Channel.Name}, 'Cz') | strcmpi({ChannelMat.Channel.Name}, 'E36'));
    iT7 = find(strcmpi({ChannelMat.Channel.Name}, 'T7') | strcmpi({ChannelMat.Channel.Name}, 'T3') | strcmpi({ChannelMat.Channel.Name}, 'E45'));
    iT8 = find(strcmpi({ChannelMat.Channel.Name}, 'T8') | strcmpi({ChannelMat.Channel.Name}, 'T4') | strcmpi({ChannelMat.Channel.Name}, 'E108'));

    if isempty(iCz) || isempty(iT7) || isempty(iT8)
        error('BST_GEODESIC_RESCUE: Nomenclature mismatch. Could not find Cz, T7, or T8.');
    end

    Cz_raw = ChannelMat.Channel(iCz(1)).Loc(:,1)';
    T7_raw = ChannelMat.Channel(iT7(1)).Loc(:,1)';
    T8_raw = ChannelMat.Channel(iT8(1)).Loc(:,1)';

    % 2. Absolute Unit Forcing (The V21.8 Hyper-Scale Empirical Shield)
    emp_rad = norm(Cz_raw - T7_raw);
    if emp_rad < 1.5, scale_f = 1.0;
    elseif emp_rad >= 1.5 && emp_rad < 30.0, scale_f = 0.01;
    elseif emp_rad >= 30.0 && emp_rad < 500.0, scale_f = 0.001;
    else, scale_f = 0.000001; end % Hyper-Scale Catch

    Cz_m = Cz_raw * scale_f; T7_m = T7_raw * scale_f; T8_m = T8_raw * scale_f;

    % 3. Resolve weights: cohort preset (if set) > UI rho/beta values
    cohort_tag = sProcess.options.cohort.Value;
    rho_val    = sProcess.options.rho.Value;
    beta_val   = sProcess.options.beta.Value;
    Dstd_val   = 0.1388;
    if ~isempty(cohort_tag)
        preset = load_cohort_preset(cohort_tag);
        rho_val = preset.rho; beta_val = preset.beta; Dstd_val = preset.D_standard;
    end

    % 4. Execute Geodesic Core Engine (Meter-Space Output)
    try
        [pLHJ_m, pRHJ_m, ~] = geodesic_rescue(Cz_m, T7_m, T8_m, ...
            'rho', rho_val, 'beta', beta_val, 'D_standard', Dstd_val);
        
        % 4. Brainstorm SCS mapping (Using Subject's structural transform)
        [sSubject, ~] = bst_get('Subject', sInputs.SubjectName);
        if isempty(sSubject.iAnatomy), error('No MRI loaded.'); end
        
        MriFile = sSubject.Anatomy(sSubject.iAnatomy).FileName;
        sMri = in_mri_bst(MriFile);
        
        % Brainstorm handles the MNI(RAS) to SCS(ALS) transformation automatically
        Pred_SCS_L = cs_convert(sMri, 'mni', 'scs', pLHJ_m);
        Pred_SCS_R = cs_convert(sMri, 'mni', 'scs', pRHJ_m);
        
        % Overwrite LPA/RPA with pristine Geodesic anchors
        sMri.SCS.LPA = Pred_SCS_L;
        sMri.SCS.RPA = Pred_SCS_R;
        
        bst_save(file_fullpath(MriFile), sMri, 'v7');
        panel_protocols('UpdateTree');
        bst_report('Info', sProcess, sInputs, sprintf('SUCCESS: Rescued origins (Scale factor: %.6f)', scale_f));
    catch ME
        bst_report('Error', sProcess, sInputs, ME.message);
    end
end
