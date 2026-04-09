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

% DIAGNOSTIC_MRI_NCS.m
% Check the NCS registration structure for Hong Kong

mri_path = 'C:\MoBI_Research\brainstorm_db\ds004718\anat\ds004718_sub-HK001\subjectimage_sub-HK001_desc-normalized_T1w.mat';
sMri = load(mri_path);

fprintf('--- MRI DIAGNOSTIC: %s ---\n', mri_path);
if isfield(sMri, 'NCS')
    disp('NCS Field found.');
    disp(sMri.NCS);
    if isfield(sMri.NCS, 'y')
        fprintf('Deformation field y size: [%s]\n', num2str(size(sMri.NCS.y)));
    else
        disp('Deformation field y is MISSING (Linear Only).');
    end
else
    disp('CRITICAL: NCS Field is MISSING.');
end

exit;
