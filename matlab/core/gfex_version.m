function v = gfex_version()
%GFEX_VERSION  Return the current GFEX-EEG toolbox version string.
%
%   v = GFEX_VERSION() returns the string identifying the installed
%   GFEX-EEG toolbox release (e.g., '1.1.5'). This is the single source
%   of truth for the MATLAB side of the toolbox — file headers and
%   wrappers refer to this function rather than hard-coding versions.
%
%   See also: gfex_version_atleast.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Version constant
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
%     See also CITATION.cff in the repository root (machine-readable).
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
%     See LICENSE in the repository root.
% ==============================================================================
    v = '1.1.5';
end
