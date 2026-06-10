function eegplugin_geodesic_rescue(fig, try_strings, catch_strings)
%EEGPLUGIN_GEODESIC_RESCUE  EEGLAB plugin loader for GFEX-EEG geodesic rescue.
%
%   Auto-discovered by EEGLAB on startup; installs a "Geodesic Origin
%   Rescue" menu item under Tools > Locate dipoles using DIPFIT (or
%   directly under Tools if DIPFIT is not available). Clicking the menu
%   item invokes pop_geodesic_rescue on the active EEG dataset to
%   generate missing LHJ / RHJ fiducials.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — EEGLAB plugin loader (MATLAB)
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

    % Define the command that executes the pop_ function
    cmd = '[EEG, LASTCOM] = pop_geodesic_rescue(EEG);';
    
    % Ensure channel locations are loaded before enabling
    cmd = [ try_strings.check_chanlocs cmd catch_strings.add_to_hist ];
    
    % Add menu item under "Tools > Locate dipoles using DIPFIT"
    % or directly under Tools if DIPFIT is not found
    toolsmenu = findobj(fig, 'tag', 'tools');
    dipfitmenu = findobj(toolsmenu, 'label', 'Locate dipoles using DIPFIT');
    
    if ~isempty(dipfitmenu)
        uimenu(dipfitmenu, 'label', 'Geodesic Origin Rescue', ...
               'callback', cmd, 'separator', 'on', 'userdata', 'startup:off;chanloc:on');
    else
        uimenu(toolsmenu, 'label', 'Geodesic Origin Rescue', ...
               'callback', cmd, 'userdata', 'startup:off;chanloc:on'); 
    end
end
