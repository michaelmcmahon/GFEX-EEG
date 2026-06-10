"""GFEX-EEG: MNE-Python black-box wrapper

Single-function entry point that appends GFEX-EEG-predicted LHJ/RHJ
coordinates to an ``mne.io.Raw`` object. Handles parity restoration
(ALS / BIDS-Brainstorm -> RAS) on the digitizer points and dispatches
to the core geodesic engine with optional cohort preset + Tier 1.5
MLP residual correction.

Usage
-----
    from geodesic_rescue_py import apply_geodesic_rescue
    raw_out = apply_geodesic_rescue(raw, cohort='LEMON_Polhemus_Adult')

Notes
-----
Coordinates assumed to be in metres (MNE standard). Weight precedence:
explicit rho/beta/D_standard kwargs > cohort preset > hard default
(LEMON_Polhemus_Adult).

Citation
--------
If you use this code in research, please cite both the software archive
and the accompanying manuscript:

    [Software]
    McMahon, M., Schukat, M., & Barrett, E. (2026).
    GFEX-EEG Toolbox [Software].
    Zenodo. https://doi.org/10.5281/zenodo.20580899

    [Paper]
    McMahon, M., Schukat, M., & Barrett, E. (Submitted).
    GFEX-EEG: Geodesic recovery of anatomical fiducials for MRI-free
    EEG source imaging.

See also CITATION.cff in the repository root (machine-readable).
Current version: see ``geodesic_rescue_py.__version__`` (single source
of truth in ``__init__.py``).

Repository
----------
https://github.com/michaelmcmahon/GFEX-EEG
Issues: https://github.com/michaelmcmahon/GFEX-EEG/issues

License
-------
MIT — see LICENSE in the repository root.
"""
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway

import mne
import numpy as np
from .core import GeodesicRescue

def apply_geodesic_rescue(raw, rho=None, beta=None, parity='RAS',
                          D_standard=None, cohort=None, verbose=True):
    """
    Appends predicted LHJ/RHJ (LPA/RPA) coordinates to mne.io.Raw.
    Coordinates are assumed to be in meters (MNE standard).
    If parity is 'ALS' or 'BIDS-Brainstorm', automatically swaps [-Y, X, Z] to restore RAS.

    Weight precedence: explicit rho/beta/D_standard kwargs >  cohort preset >
    hard default (LEMON_Polhemus_Adult, ~16.8 mm held-out on LEMON HTJ GT).

    Parameters
    ----------
    cohort : str, optional
        Cohort tag from weight_zoo.json (e.g. 'LEMON_Polhemus_Adult',
        'EGI_HydroCel_256_Adult', 'lemon', 'default'). See
        ``load_cohort_preset`` for the full list.
    """
    if not isinstance(raw, mne.io.BaseRaw):
        raise ValueError("Input must be an MNE Raw object.")
    
    # 1. Apply global parity swap to ALL digitizer points if ALS
    if parity.upper() in ['ALS', 'BIDS-BRAINSTORM']:
        with raw.info._unlock():
            for d in raw.info['dig']:
                # ALS to RAS: [-Y, X, Z]
                d['r'] = np.array([-d['r'][1], d['r'][0], d['r'][2]])
    
    # 2. Identify Anchors
    def find_anchor(names):
        for d in raw.info['dig']:
            # Check ch_names or labels
            # In MNE, fiducials are usually FIFFV_POINT_LPA, etc.
            # But here we look for EEG channels as anchors
            pass
        
        # Look in channel names
        ch_names = raw.ch_names
        for name in names:
            for ch in ch_names:
                # Remove common prefixes
                clean_ch = ch.upper().replace('EEG', '').replace('E', '')
                if clean_ch == name.upper() or ch.upper() == name.upper():
                    idx = ch_names.index(ch)
                    # Get position from info['chs']
                    pos = raw.info['chs'][idx]['loc'][:3]
                    if np.all(pos == 0):
                        continue
                    return pos
        return None

    Cz_m = find_anchor(['Cz', '80', '18', '36'])
    T7_m = find_anchor(['T7', 'T3', '45', '13', '44', '46'])
    T8_m = find_anchor(['T8', 'T4', '108', '14', '77', '102'])

    if Cz_m is None or T7_m is None or T8_m is None:
        raise RuntimeError("Could not find Cz, T7, and T8 anchors in the Raw object.")

    if parity.upper() in ['ALS', 'BIDS-BRAINSTORM']:
        Cz_m = np.array([-Cz_m[1], Cz_m[0], Cz_m[2]])
        T7_m = np.array([-T7_m[1], T7_m[0], T7_m[2]])
        T8_m = np.array([-T8_m[1], T8_m[0], T8_m[2]])

    if verbose:
        print(f"Anchors found: Cz={Cz_m}, T7={T7_m}, T8={T8_m}")

    # 2. Execute Rescue (weight resolution happens in GeodesicRescue.__init__)
    rescuer = GeodesicRescue(rho=rho, beta=beta, D_standard=D_standard, cohort=cohort)
    pLHJ, pRHJ = rescuer.rescue(Cz_m, T7_m, T8_m)

    if verbose:
        print(f"Predicted LHJ (LPA): {pLHJ}")
        print(f"Predicted RHJ (RPA): {pRHJ}")

    # 3. Append to info['dig']
    # MNE Fiducial IDs:
    # FIFFV_POINT_CARDINAL = 1
    # Cardinal IDs: LPA=1, NASION=2, RPA=3
    
    def update_fiducial(ident, pos, name):
        found = False
        for i, d in enumerate(raw.info['dig']):
            if d['kind'] == 1 and d['ident'] == ident:
                raw.info['dig'][i]['r'] = pos.astype(np.float32)
                found = True
                if verbose: print(f"Updated existing {name} in dig structure.")
                break
        if not found:
            new_dig = {
                'kind': 1,
                'ident': ident,
                'r': pos.astype(np.float32),
                'coord_frame': 4 # FIFFV_COORD_HEAD
            }
            raw.info.set_montage(mne.channels.make_dig_montage(
                lpa=pLHJ if ident==1 else None,
                rpa=pRHJ if ident==3 else None,
                coord_frame='head'
            ), on_missing='ignore')
            # Actually, manually appending is safer if we want to preserve existing dig
            raw.info['dig'].append({
                'kind': 1,
                'ident': ident,
                'r': pos.astype(np.float32),
                'coord_frame': 4
            })
            if verbose: print(f"Appended new {name} to dig structure.")

    update_fiducial(1, pLHJ, 'LPA')
    update_fiducial(3, pRHJ, 'RPA')

    return raw
