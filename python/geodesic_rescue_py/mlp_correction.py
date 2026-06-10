"""GFEX-EEG: Tier 1.5 residual-correction MLP (Python)

Residual-learning MLP correction on top of the geodesic prediction.
Single-hidden-layer network (15-input -> 64 ReLU -> 6-output, 1,414
parameters) trained against per-subject MRI-tagged HTJ ground truth.

Architecture
------------
    Inputs : (Cz, T7, T8, pL_geo, pR_geo) in RAS metres
    Outputs: (dL, dR) correction vectors in metres
    Final  : pL_corrected = pL_geo + dL,  pR_corrected = pR_geo + dR

Inference is pure matrix multiplication — no ML-framework dependency at
runtime (NumPy only). Weights are cached per-path to avoid repeated
disk IO.

Three production cohort presets ship: ``LEMON_Polhemus_Adult`` (Polhemus,
4.76 mm held-out), ``WH_Neuromag70_Adult`` (Neuromag, 8.16 mm held-out),
``CapTrak_Adult`` (BrainProducts CapTrak, 8.52 mm leave-one-out).

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

import os
import numpy as np
import scipy.io as sio

_CACHE = {}


def _load_mlp(mlp_path):
    mlp_path = os.path.abspath(mlp_path)
    if mlp_path in _CACHE:
        return _CACHE[mlp_path]
    if not os.path.isfile(mlp_path):
        raise FileNotFoundError(f'MLP weights file not found: {mlp_path}')
    raw = sio.loadmat(mlp_path)
    cache = {
        'W1':     np.asarray(raw['W1'],     dtype=float),
        'b1':     np.asarray(raw['b1'],     dtype=float).flatten(),
        'W2':     np.asarray(raw['W2'],     dtype=float),
        'b2':     np.asarray(raw['b2'],     dtype=float).flatten(),
        'x_mean': np.asarray(raw['x_mean'], dtype=float).flatten(),
        'x_std':  np.asarray(raw['x_std'],  dtype=float).flatten(),
        'y_mean': np.asarray(raw['y_mean'], dtype=float).flatten(),
        'y_std':  np.asarray(raw['y_std'],  dtype=float).flatten(),
    }
    _CACHE[mlp_path] = cache
    return cache


def apply_mlp_correction(Cz, T7, T8, pL_geo, pR_geo, mlp_path):
    """Return (pL_corrected, pR_corrected, info_dict) in metres.

    All inputs are length-3 vectors in RAS metres. info dict contains
    'dL', 'dR', 'dL_magnitude_mm', 'dR_magnitude_mm'.
    """
    mlp = _load_mlp(mlp_path)

    # 1. Build input feature vector (1D, 15)
    x = np.concatenate([
        np.asarray(Cz).flatten(),
        np.asarray(T7).flatten(),
        np.asarray(T8).flatten(),
        np.asarray(pL_geo).flatten(),
        np.asarray(pR_geo).flatten(),
    ])

    # 2. Normalise input
    x_norm = (x - mlp['x_mean']) / mlp['x_std']

    # 3. Hidden layer (ReLU)
    h = x_norm @ mlp['W1'] + mlp['b1']
    h = np.maximum(h, 0.0)

    # 4. Output layer (linear)
    y_norm = h @ mlp['W2'] + mlp['b2']

    # 5. Denormalise output
    y = y_norm * mlp['y_std'] + mlp['y_mean']
    dL = y[0:3]
    dR = y[3:6]

    pL_corr = np.asarray(pL_geo).flatten() + dL
    pR_corr = np.asarray(pR_geo).flatten() + dR
    info = {
        'dL': dL,
        'dR': dR,
        'dL_magnitude_mm': float(np.linalg.norm(dL) * 1000),
        'dR_magnitude_mm': float(np.linalg.norm(dR) * 1000),
    }
    return pL_corr, pR_corr, info
