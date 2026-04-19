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

function varargout = process_geodesic_rescue( varargin )
% PROCESS_GEODESIC_RESCUE (V21.8 - Invincible Brainstorm Wrapper)
% Includes: Hyper-Scale Empirical Shield and Brainstorm SCS mapping.

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

    % 3. Execute Geodesic Core Engine (Meter-Space Output)
    try
        [pLHJ_m, pRHJ_m, ~] = geodesic_rescue(Cz_m, T7_m, T8_m, 'rho', sProcess.options.rho.Value, 'beta', sProcess.options.beta.Value);
        
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
