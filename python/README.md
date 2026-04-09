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

# Geodesic Rescue Python (geodesic-rescue-py)

Python implementation of the Geodesic Rescue Toolbox for EEG fiducial extrapolation.

## Installation

```bash
pip install .
```

## Usage with MNE-Python

```python
import mne
from geodesic_rescue_py import apply_geodesic_rescue

# Load your data
raw = mne.io.read_raw_fif("my_data.fif")

# Apply rescue
apply_geodesic_rescue(raw)

# Now raw.info['dig'] contains predicted LPA/RPA
```
