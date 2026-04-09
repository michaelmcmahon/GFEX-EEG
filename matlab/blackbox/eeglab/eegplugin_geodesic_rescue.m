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

function eegplugin_geodesic_rescue(fig, try_strings, catch_strings)
% eegplugin_geodesic_rescue - EEGLAB plugin to generate missing LHJ/RHJ fiducials

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
