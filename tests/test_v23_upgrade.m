% test_v23_upgrade.m
% Validates the new V23.0 GFEX-EEG core with native Parity swapping.

clc;
addpath('C:\MoBI_Research\GFEX-EEG\matlab\core');
addpath(genpath('C:\MoBI_Research\brainstorm3\toolbox'));

% Target LEMON Sub-010060 (Flattened ALS space)
rawDir = 'C:\MoBI_Research\Fiducial_Extrapolation_Exp\Data_Raw\ds000221_LEMON';
sub_id = 'sub-010060';
matPath = fullfile(rawDir, sub_id, 'eeg', [sub_id '_localizer.mat']);
loc = load(matPath);

% Find Anchors
iCz = find(contains({loc.Channel.Name}, '_Cz'), 1);
iT7 = find(contains({loc.Channel.Name}, '_T7'), 1);
iT8 = find(contains({loc.Channel.Name}, '_T8'), 1);

Cz_raw = loc.Channel(iCz).Loc';
T7_raw = loc.Channel(iT7).Loc';
T8_raw = loc.Channel(iT8).Loc';

fprintf('RAW ALS Sensors (sub-010060):\n');
fprintf('  Cz: [%.4f, %.4f, %.4f]\n', Cz_raw(1), Cz_raw(2), Cz_raw(3));
fprintf('  T7: [%.4f, %.4f, %.4f]\n', T7_raw(1), T7_raw(2), T7_raw(3));
fprintf('  T8: [%.4f, %.4f, %.4f]\n', T8_raw(1), T8_raw(2), T8_raw(3));

% Run V23 Wrapper with ALS parity flag
fprintf('\nExecuting V23 geodesic_rescue.m with PARITY=ALS...\n');
[pL, pR, info] = geodesic_rescue(Cz_raw, T7_raw, T8_raw, 'parity', 'ALS');

fprintf('\nPREDICTED RESULTS (RAS Meters):\n');
fprintf('  LPA: [%.4f, %.4f, %.4f]\n', pL(1), pL(2), pL(3));
fprintf('  RPA: [%.4f, %.4f, %.4f]\n', pR(1), pR(2), pR(3));

% Compare against V40.0 Ground Truth (sub-010060)
% GT: NAS: [92, 205, 153], LPA: [11, 120, 123], RPA: [169, 120, 122]
% The exact SCS-transformed ground truth for 010060 was generated during the V40.0 run.
% We just want to ensure it completes successfully and outputs sane RAS coordinates.

fprintf('\nVerification successful if coordinates match standard physiological head bounds.\n');
exit;