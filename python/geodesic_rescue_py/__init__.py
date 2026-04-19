# /*******************************************************************************
# * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
# * Version: 1.0.0
# * Repository: https://github.com/michaelmcmahon/GFEX-EEG
# * License:  MIT License
# * Authors: Michael McMahon / University of Galway
# * DOI: [If available]
# * Date: 2026
# *
# * [License Text or link to License file]
# *******************************************************************************/

from .core import GeodesicRescue
from .mne_wrapper import apply_geodesic_rescue
from .tuner import GeodesicTuner
from .weight_zoo import load_cohort_preset
from .mlp_correction import apply_mlp_correction

__version__ = "1.1.0.dev0"
