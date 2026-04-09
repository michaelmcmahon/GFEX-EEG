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

# GFEX-EEG Toolbox (V1.0.0)

**A Hardware-Agnostic Spatial Extrapolation Engine for EEG Fiducials**

GFEX-EEG solves the "MRI-Free, Fiducial-Free" paradox for high-density EEG systems (like EGI and CapTrak) and recovers floating origins from legacy datasets (like Wakeman-Henson). It uses a Far & Near Optimization (FNO) Metaheuristic paired with a Procrustes-Dijkstra Manifold Engine to extrapolate exact spatial locations for the Left and Right Pre-Auricular Points (LPA/RPA) using only the structural geometry of the Cz, T7, and T8 electrodes.

## Distribution Architecture

This toolbox is distributed as a single integrated package logically separated into three tiers. This design strictly adheres to the DRY (Don't Repeat Yourself) principle, ensuring the "Black-Box" ease-of-use wrappers and the Advanced "Tuning" mode both query the exact same physical engine.

```text
GFEX-EEG/
├── data/
│   └── ICBM152_scalp.mat          <- Shared Euclidean manifold substrate
├── matlab/
│   ├── core/                      <- Shared Physics Engine
│   │   ├── geodesic_rescue.m
│   │   └── predict_helix_tragus_junctions_fno.m
│   ├── blackbox/                  <- Standard "Black-Box" Wrappers
│   │   ├── brainstorm/
│   │   ├── eeglab/
│   │   └── fieldtrip/
│   └── tuning/                    <- Advanced Metaheuristic Mode
│       ├── geodesic_fno_tuner.m
│       └── verify_geodesic_tuner.m
├── python/                        <- Python Native Library (geodesic-rescue-py)
│   ├── geodesic_rescue_py/
│   │   ├── core.py                <- Shared physics (Scipy/Numpy)
│   │   ├── mne_wrapper.py         <- MNE-Python "Black-Box"
│   │   └── tuner.py               <- GeodesicTuner Class
│   └── setup.py
├── reproducibility/               <- Experimental Code Pipeline (Paper Validation)
└── README.md
```

## Tier 1: "Black-Box" Wrappers (For Standard Adult Cohorts)

For 95% of use cases, researchers can rely on the pre-trained weights (`rho = 0.0054`, `beta = 1.1885`). The wrappers automatically apply the **Hyper-Scale Catch** to neutralize BIDS double-scaling traps and execute **RAS-to-ALS** axis transposition.

### EEGLAB
Simply add `matlab/blackbox/eeglab` to your MATLAB path. The plugin will appear under `Tools > Geodesic Origin Rescue`.
Or call it programmatically:
```matlab
[EEG, ~] = pop_geodesic_rescue(EEG);
```

### FieldTrip
```matlab
cfg = [];
elec_rescued = ft_geodesic_rescue(cfg, elec_raw);
```

### MNE-Python
```bash
cd python
pip install .
```
```python
import mne
from geodesic_rescue_py import apply_geodesic_rescue

raw = mne.io.read_raw_fif("my_data.fif")
raw_rescued = apply_geodesic_rescue(raw)
```

## Tier 2: Metaheuristic Tuning Mode (For Advanced Users)

While the pre-trained weights cover the vast majority of adult human heads, researchers working with highly specialized hardware (e.g., neonatal caps) or pediatric cohorts with vastly different cranial proportions will need to recalculate the optimal parameters.

If you possess a small "training subset" (e.g., 5 to 10 subjects) where you actually have structural MRIs and have manually tagged the true LHJ/RHJ, you can use the dual-phase Far and Near Optimization (FNO) engine to derive bespoke bounds.

### MATLAB Tuner
```matlab
% training_data requires: .Cz, .T7, .T8, .gtLPA, .gtRPA (in Meters)
[optimal_rho, optimal_beta, info] = geodesic_fno_tuner(training_data);
```

### Python Tuner
```python
from geodesic_rescue_py import GeodesicTuner

tuner = GeodesicTuner(training_data)
results = tuner.tune()
print(f"Optimal Rho: {results['optimal_rho']}")
```

Once derived, you can pass these optimal parameters into the Black-Box wrappers to rescue the remaining 90% of your dataset that lacks MRIs.

---
*Authors: Michael McMahon / University of Galway*
