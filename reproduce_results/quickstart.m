% =============================================================================
% GFEX-EEG quickstart (MATLAB)
% =============================================================================
% Sanity-check the toolbox end-to-end without downloading any external dataset.
%
% Uses the shipped template Cz/T7/T8 anchors (canonical RAS-metre positions on
% the ICBM152 scalp manifold) as a synthetic subject input, then exercises:
%
%   1. Pure geodesic prediction        (Tier 1)
%   2. LEMON_Polhemus_Adult MLP        (Tier 1.5)
%   3. WH_Neuromag70_Adult MLP         (Tier 1.5)
%   4. CapTrak_Adult MLP               (Tier 1.5)
%
% For each mode the predicted Helix-Tragus Junction coordinates are compared
% against hardcoded golden values (machine-precision tolerance: 1e-5 m =
% 0.01 mm). All four assertions should PASS on a clean install.
%
% This script demonstrates that the toolbox loads, runs, and produces
% deterministic bit-identical output across runs. It does NOT reproduce the
% per-cohort accuracy numbers (4.76 / 8.16 / 8.52 mm) reported in the paper:
% those additionally require downloading the public LEMON / Wakeman-Henson /
% ds004024 datasets per the manuscript's Data Availability statement.
%
% Usage (from MATLAB, with the toolbox on the path):
%   >> cd <toolbox-root>/reproduce_results
%   >> quickstart
%
% Expected runtime: < 1 second.
% =============================================================================

% Add toolbox core to MATLAB path
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'matlab', 'core'));

% ---- Synthetic subject input: shipped template Cz/T7/T8 anchors ------------
Cz = [0.011240, 0.025921, 0.141134];
T7 = [-0.089174, -0.001327, -0.006348];
T8 = [0.096880, -0.014286, -0.005819];

% ---- Golden expected outputs (must match expected_output.json) -------------
EXPECTED = struct();
EXPECTED.pure_geodesic.LHJ = [-0.09101699973067885, -0.03453133203031224, -0.05497997557054298];
EXPECTED.pure_geodesic.RHJ = [ 0.0901168635406814,  -0.015252983816377925, -0.062139438734905815];
EXPECTED.lemon_mlp.LHJ     = [-0.07709191538808385, -0.024460965635807783, -0.0014766207622546393];
EXPECTED.lemon_mlp.RHJ     = [ 0.0850239015021078,  -0.00853983954755081,  -0.004513472947956558];
EXPECTED.wh_mlp.LHJ        = [-0.12843714410646706, -0.05771261485109747,  -0.028998996630177232];
EXPECTED.wh_mlp.RHJ        = [ 0.1148258931103476,  -0.022229687867550817, -0.03979050417794089];
EXPECTED.captrak_mlp.LHJ   = [-0.10381647373226766, -0.0232966729954635,   -0.02303326264563576];
EXPECTED.captrak_mlp.RHJ   = [ 0.1099652541965088,  -0.0063765445063921275, -0.02771504842231108];

TOL = 1e-5;  % 0.01 mm

fprintf('GFEX-EEG quickstart (MATLAB) -- end-to-end sanity check\n');
fprintf('========================================================================\n\n');
fprintf('Input anchors (RAS metres):\n');
fprintf('  Cz: [%+.6f, %+.6f, %+.6f]\n', Cz(1), Cz(2), Cz(3));
fprintf('  T7: [%+.6f, %+.6f, %+.6f]\n', T7(1), T7(2), T7(3));
fprintf('  T8: [%+.6f, %+.6f, %+.6f]\n', T8(1), T8(2), T8(3));
fprintf('\n');

% ---- Demo 1: pure geodesic --------------------------------------------------
fprintf('Demo 1 -- Pure geodesic (Tier 1)\n');
[lhj, rhj] = geodesic_rescue(Cz, T7, T8, 'cohort', 'LEMON_Polhemus_Adult');
check_('LEMON cohort, pure geodesic', lhj, rhj, EXPECTED.pure_geodesic, TOL);
fprintf('\n');

% ---- Demo 2: LEMON MLP ------------------------------------------------------
fprintf('Demo 2 -- LEMON_Polhemus_Adult MLP correction (Tier 1.5)\n');
[lhj, rhj] = geodesic_rescue(Cz, T7, T8, ...
    'cohort', 'LEMON_Polhemus_Adult', 'mlp_correction', true);
check_('LEMON_Polhemus_Adult + MLP', lhj, rhj, EXPECTED.lemon_mlp, TOL);
fprintf('\n');

% ---- Demo 3: WH MLP ---------------------------------------------------------
fprintf('Demo 3 -- WH_Neuromag70_Adult MLP correction (Tier 1.5)\n');
[lhj, rhj] = geodesic_rescue(Cz, T7, T8, ...
    'cohort', 'WH_Neuromag70_Adult', 'mlp_correction', true);
check_('WH_Neuromag70_Adult + MLP', lhj, rhj, EXPECTED.wh_mlp, TOL);
fprintf('\n');

% ---- Demo 4: CapTrak MLP ----------------------------------------------------
fprintf('Demo 4 -- CapTrak_Adult MLP correction (Tier 1.5)\n');
[lhj, rhj] = geodesic_rescue(Cz, T7, T8, ...
    'cohort', 'CapTrak_Adult', 'mlp_correction', true);
check_('CapTrak_Adult + MLP', lhj, rhj, EXPECTED.captrak_mlp, TOL);
fprintf('\n');

fprintf('========================================================================\n');
fprintf('All 4 demos PASSED. Toolbox is functioning end-to-end.\n\n');
fprintf('Note: this verifies deterministic engine behaviour against golden\n');
fprintf('values. Reproducing the per-cohort accuracy numbers (4.76 mm LEMON,\n');
fprintf('8.16 mm Wakeman-Henson, 8.52 mm CapTrak) reported in the manuscript\n');
fprintf('additionally requires downloading the public datasets named in the\n');
fprintf('paper''s Data Availability statement.\n');


% =============================================================================
function check_(label, lhj, rhj, expected, tol)
    delta_l = norm(lhj - expected.LHJ);
    delta_r = norm(rhj - expected.RHJ);
    pass_l = delta_l < tol;
    pass_r = delta_r < tol;
    if pass_l && pass_r
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('  %-30s  delta L = %.6f mm, delta R = %.6f mm   [%s]\n', ...
        label, delta_l*1000, delta_r*1000, status);
    fprintf('    LHJ predicted: [%+.6f, %+.6f, %+.6f] m\n', lhj(1), lhj(2), lhj(3));
    fprintf('    RHJ predicted: [%+.6f, %+.6f, %+.6f] m\n', rhj(1), rhj(2), rhj(3));
    assert(pass_l && pass_r, '%s failed golden-value check', label);
end
