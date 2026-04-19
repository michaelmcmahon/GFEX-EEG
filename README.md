# GFEX-EEG Toolbox (V1.0.0)

**A Hardware-Agnostic Spatial Extrapolation Engine for EEG Fiducials**

GFEX-EEG solves the "MRI-Free, Fiducial-Free" paradox for high-density EEG systems (like EGI and CapTrak) and recovers floating origins from legacy datasets (like Wakeman-Henson). It uses a Far & Near Optimization (FNO) Metaheuristic paired with a Procrustes-Dijkstra Manifold Engine to extrapolate **Helix-Tragus Junction (HTJ)** coordinates — a geometrically unambiguous anatomical landmark — using only the positions of the Cz, T7, and T8 electrodes. For downstream compatibility, predicted HTJ coordinates are written into the standard `LPA`/`RPA` slots of MNE-Python Raw / EEGLAB / FieldTrip / Brainstorm data structures; existing source-imaging pipelines consume them as anatomical fiducials without modification.

## Distribution Architecture

This toolbox is distributed as a single integrated package logically separated into three tiers. This design strictly adheres to the DRY (Don't Repeat Yourself) principle, ensuring the "Black-Box" ease-of-use wrappers and the Advanced "Tuning" mode both query the exact same physical engine.

```text
GFEX-EEG/
├── data/
│   └── ICBM152_scalp.mat          <- Shared Euclidean manifold substrate
├── matlab/
│   ├── core/                      <- Shared Physics Engine
│   ├── blackbox/                  <- Standard "Black-Box" Wrappers (EEGLAB, FieldTrip, Brainstorm)
│   └── tuning/                    <- Advanced Metaheuristic Mode (FNO Tuner)
├── python/                        <- Python Native Library (geodesic-rescue-py)
│   └── geodesic_rescue_py/        <- MNE-Python Wrapper & Tuner
├── reproducibility/               <- Experimental Code Pipeline (Paper Validation)
└── README.md
```

## Forensic Data Warning: The "Brainstorm Z=0 Flatline" and "Raw Voxel Trap"
When extracting ground truth from legacy neuroinformatics pipelines, **always ensure your data hasn't been structurally flattened or grid-locked.** 

During our analysis of the LEMON dataset, we proved that standard Brainstorm BIDS exports mathematically coerce the 3D Polhemus coordinates to the `Z=0` plane, effectively **erasing the subject's physiological head tilt**. Furthermore, default Brainstorm MRI coordinates are frequently stored as **unscaled, un-tilted Voxel Indices**. 

**The Fix:** Always extract True World Space coordinates (`cs_convert` to 'world' in Brainstorm) to restore the subject's physiological tilt and 1.0mm scaling before running any FNO Tuner derivations. Failure to do so will result in an Orthogonal Parity Collision.

## Tier 1: "Black-Box" Wrappers (For Standard Adult Cohorts)

For 95% of use cases, researchers can rely on the LEMON-tuned weights (`rho = 0.248383`, `beta = 0.235926`, tuned against manually-annotated HTJ ground truth on N=10 LEMON subjects, held-out mean error 16.82 mm on N=10 independent subjects). The wrappers automatically apply the **Hyper-Scale Catch** to neutralize BIDS double-scaling traps and execute **RAS-to-ALS** axis transposition.

For non-standard hardware (EGI, dense pediatric caps, custom Polhemus protocols) researchers are advised to retune via the FNO metaheuristic on a small cohort-specific training set — see Tier 2 below.

### Validation & Accuracy

On a held-out set of N=10 LEMON subjects, the Tier 1 weights achieve **16.82 mm mean error (SD 4.03 mm)** against manually-annotated HTJ ground truth.

Two caveats researchers should know when interpreting or reproducing public-dataset benchmarks:

- **HAD (ds007353) and NOD (ds005811-ds005810) ship "packed" `*_electrodes.tsv` files.** The same Cz/T7/T8 template is replicated across all subjects and sessions (82 of 83 files share a single SHA1 hash in our audit). Per-subject accuracy numbers on these datasets therefore reflect anatomical scatter of the MRI HTJ landmark at a *single fixed algorithm output*, not genuine per-subject generalization. Treat as one effective EEG input point, not N=46.

- **Cohort-specific FNO retune recommended** when deploying to non-Polhemus digitization (CapTrak, EGI). The Tier 1 LEMON-tuned weights transfer approximately; per-cohort retune via Tier 2 typically improves accuracy by several millimetres.

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

Researchers working with highly specialized hardware (e.g., neonatal caps) or pediatric cohorts will need to recalculate the optimal parameters. Using a small "training subset" (e.g., 5 subjects) with known MRIs, the FNO engine can derive bespoke bounds.

### MATLAB Tuner
```matlab
[optimal_rho, optimal_beta, info] = geodesic_fno_tuner(training_data);
```

### Python Tuner
```python
from geodesic_rescue_py import GeodesicTuner
tuner = GeodesicTuner(training_data)
results = tuner.tune()
```

---
**Authors:** Michael McMahon / University of Galway  
**Repository:** [https://github.com/michaelmcmahon/GFEX-EEG](https://github.com/michaelmcmahon/GFEX-EEG)  
**License:** MIT License  
**Date:** 2026
