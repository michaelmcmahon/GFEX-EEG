# /*******************************************************************************
# * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
# * Version: 1.0.0
# * Repository: https://github.com/michaelmcmahon/GFEX-EEG
# * License:  MIT License
# * Authors: Michael McMahon / University of Galway
# *******************************************************************************/

"""ICBM152 template positions for standard 10-20 electrodes.

All positions are mesh-derived from ICBM152_scalp.mat via 10-20 geodesic
arc proportions (Jasper 1958 / Chatrian 1985). Computation script:
    Fiducial_Extrapolation_Exp/Scripts/c_Compute_Template_Anchors.m
Derivation log:
    Fiducial_Extrapolation_Exp/Results/c_Compute_Template_Anchors_20260419.log

Method: Nz/Iz located geometrically (anterior/posterior scalp-base on the
X=0 midline at ear-level Z). Sagittal arc Nz->Cz->Iz traced via Dijkstra on
scalp-mesh edges. Midline anchors placed at canonical 10-20 proportions
(Fpz=10%, Fz=30%, Pz=70%, Oz=90%). Lateral anchors placed along Fpz->T7
and Fpz->T8 oblique arcs at 50% (F7/F8) and 75% (FT7/FT8).

Cz, T7, T8 retain the originally shipped engine reference positions.
"""

import numpy as np

TEMPLATE_ANCHORS = {
    'Cz':  ( 0.011240,  0.025921,  0.141134),
    'T7':  (-0.089174, -0.001327, -0.006348),
    'T8':  ( 0.096880, -0.014286, -0.005819),
    'Nz':  (-0.012142,  0.087680,  0.010832),
    'Iz':  (-0.009680, -0.088407,  0.009500),
    'Fpz': (-0.010372,  0.084290,  0.048318),
    'Fz':  (-0.002715,  0.066019,  0.115912),
    'Pz':  (-0.000992, -0.068584,  0.114015),
    'Oz':  (-0.008064, -0.085547,  0.046371),
    'F7':  (-0.060655,  0.057793, -0.006690),
    'F8':  ( 0.073017,  0.058577, -0.019214),
    'FT7': (-0.080462,  0.029685, -0.020110),
    'FT8': ( 0.092837,  0.023594, -0.023749),
}


def template_anchor(label):
    """Return a 3-vector (metres, RAS) for a supported 10-20 label.

    Parameters
    ----------
    label : str
        Electrode label (case-sensitive), e.g. 'Fpz', 'Oz', 'FT7'.

    Returns
    -------
    numpy.ndarray
        Shape (3,) — (x, y, z) coordinate in metres, RAS frame.

    Raises
    ------
    KeyError
        If label is not in the supported template set.
    """
    if label not in TEMPLATE_ANCHORS:
        raise KeyError(
            f'No template position for anchor "{label}". '
            f'Supported labels: {", ".join(TEMPLATE_ANCHORS.keys())}'
        )
    return np.array(TEMPLATE_ANCHORS[label])


def supported_labels():
    """Return list of all supported anchor labels."""
    return list(TEMPLATE_ANCHORS.keys())
