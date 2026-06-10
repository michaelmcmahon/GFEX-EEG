function tf = gfex_version_atleast(required)
%GFEX_VERSION_ATLEAST  Test installed GFEX-EEG version against a requirement.
%
%   tf = GFEX_VERSION_ATLEAST(REQUIRED) returns true if the installed
%   GFEX-EEG version (from gfex_version()) is greater than or equal to
%   REQUIRED. Comparison is numeric per dotted segment — '1.1.10' is
%   correctly treated as greater than '1.1.2' (not lexically).
%
%   USAGE
%       if ~gfex_version_atleast('1.1.5')
%           error('Update GFEX-EEG to v1.1.5 or later');
%       end
%
%   See also: gfex_version.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Version comparison helper
% ------------------------------------------------------------------------------
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
% ==============================================================================
    cur = sscanf(gfex_version(), '%d.%d.%d');
    req = sscanf(required, '%d.%d.%d');
    for i = 1:min(numel(cur), numel(req))
        if cur(i) > req(i), tf = true;  return; end
        if cur(i) < req(i), tf = false; return; end
    end
    tf = numel(cur) >= numel(req);
end
