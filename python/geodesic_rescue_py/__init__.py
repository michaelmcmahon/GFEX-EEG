"""GFEX-EEG: Python package entry point

Public Python interface for the GFEX-EEG toolbox — geodesic fiducial
extrapolation for MRI-free EEG source imaging. Exposes the deterministic
geodesic walk engine, cohort-specific MLP residual correction, the
Far/Near hyperparameter tuner, and an MNE-Python-compatible wrapper.

Usage
-----
    from geodesic_rescue_py import apply_geodesic_rescue, GeodesicRescue

    # Functional one-shot (default LEMON-tuned weights)
    pL, pR = apply_geodesic_rescue(Cz, T7, T8)

    # Object form (cohort preset + opt-in MLP residual correction)
    rescuer = GeodesicRescue(cohort='LEMON_Polhemus_Adult', mlp_correction=True)
    pL, pR = rescuer.rescue(Cz, T7, T8)

Programmatic version checks
---------------------------
    from packaging.version import Version
    import geodesic_rescue_py as gfex
    if Version(gfex.__version__) < Version('1.1.5'):
        raise RuntimeError('Update GFEX-EEG to 1.1.5 or later')

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

# --- Package metadata (single source of truth for the Python side) ------------
__version__ = "1.1.6"
__author__ = "Michael McMahon, Michael Schukat, Enda Barrett"
__license__ = "MIT"
__citation_doi__ = "10.5281/zenodo.20580898"
__citation_url__ = "https://doi.org/10.5281/zenodo.20580898"

# --- Public API re-exports ---------------------------------------------------
from .core import GeodesicRescue
from .mne_wrapper import apply_geodesic_rescue
from .tuner import GeodesicTuner
from .weight_zoo import load_cohort_preset
from .mlp_correction import apply_mlp_correction
