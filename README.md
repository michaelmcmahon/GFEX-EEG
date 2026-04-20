# GFEX-EEG Toolbox (V1.1.0)

**A Hardware-Agnostic Spatial Extrapolation Engine for EEG Fiducials**

GFEX-EEG solves the "MRI-Free, Fiducial-Free" paradox for high-density EEG systems (like EGI and CapTrak) and recovers floating origins from legacy datasets (like Wakeman-Henson). It uses a Far & Near Optimization (FNO) Metaheuristic paired with a Procrustes-Dijkstra Manifold Engine to extrapolate **Helix-Tragus Junction (HTJ)** coordinates — a geometrically unambiguous anatomical landmark — using only the positions of the Cz, T7, and T8 electrodes. For downstream compatibility, predicted HTJ coordinates are written into the standard `LPA`/`RPA` slots of MNE-Python Raw / EEGLAB / FieldTrip / Brainstorm data structures; existing source-imaging pipelines consume them as anatomical fiducials without modification.

**V1.1.0 adds an optional Tier 1.5 residual-correction MLP layered on top of the geodesic prediction**, reducing LEMON held-out error from 16.82 mm to 4.76 mm (71.7% reduction). The MLP is opt-in via a single kwarg; the default pure-geodesic path is preserved bit-identically and remains tagged at `v1.0.0-pure-geodesic` for reproducibility of pre-V1.1 results. See the Tier 1.5 section below.

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

For non-standard hardware (EGI, dense pediatric caps, custom Polhemus protocols) researchers are advised to retune via the FNO metaheuristic on a small cohort-specific training set — see Tier 2 below. Contributed cohort-specific weights can be consumed via the **cohort-preset** mechanism (see next section).

### Cohort Presets (Weight Zoo)

Production weights are versioned in `data/weight_zoo.json`. Any wrapper accepts a `cohort` tag to load a preset without specifying ρ/β manually:

```matlab
% MATLAB — cohort preset (zoo lookup)
[pL, pR] = geodesic_rescue(Cz, T7, T8, 'cohort', 'LEMON_Polhemus_Adult');
% or with a named alias:
[pL, pR] = geodesic_rescue(Cz, T7, T8, 'cohort', 'lemon');

% Explicit rho/beta override the preset:
[pL, pR] = geodesic_rescue(Cz, T7, T8, 'cohort', 'LEMON_Polhemus_Adult', 'rho', 0.25);
```

```python
# Python
from geodesic_rescue_py import apply_geodesic_rescue, load_cohort_preset

raw = apply_geodesic_rescue(raw, cohort='LEMON_Polhemus_Adult')

# Inspect a preset directly:
preset = load_cohort_preset('lemon')   # {'rho': 0.248383, 'beta': 0.235926, ...}
```

**Currently shipped production entries:** `LEMON_Polhemus_Adult` (default, 16.82 mm held-out).

**Pending community contribution:** `CapTrak_Adult`, `EGI_HydroCel_256_Adult`, `Neuroscan_SynAmps_Adult` — entries are reserved in the zoo with `status: pending_tune`. Researchers with MRI-bearing cohorts on those setups can run the Tier 2 FNO tuner and submit a PR adding a production entry. See `data/weight_zoo.json` `contribution` block for the submission checklist.

### Validation & Accuracy

On a held-out set of N=10 LEMON subjects with manually-annotated HTJ ground truth:

| Mode | Mean error | SD | Improvement vs pure geodesic |
|---|---|---|---|
| Tier 1 pure geodesic | **16.82 mm** | 4.03 | — |
| Tier 1.5 (+ MLP residual correction) | **4.76 mm** | 2.06 | **71.7% reduction** |

The Tier 1.5 result survives a 4-test data-leak diligence battery (subject separation; label-shuffle permutation test with 3× degradation on scrambled labels; alternative-holdout rotation at 5.29 ± 0.62 mm across 5 random 10-subject splits; FNO-seen vs unseen within training at 0.15 mm gap). Full report: `Fiducial_Extrapolation_Exp/Results/c_MLP_Diligence_20260420.json`.

Two caveats researchers should know when interpreting or reproducing public-dataset benchmarks:

- **HAD (ds007353) and NOD (ds005811-ds005810) ship "packed" `*_electrodes.tsv` files.** The same Cz/T7/T8 template is replicated across all subjects and sessions (82 of 83 files share a single SHA1 hash in our audit). Per-subject accuracy numbers on these datasets therefore reflect anatomical scatter of the MRI HTJ landmark at a *single fixed algorithm output*, not genuine per-subject generalization. Treat as one effective EEG input point, not N=46.

- **Cohort-specific FNO retune recommended** when deploying to non-Polhemus digitization (CapTrak, EGI). The Tier 1 LEMON-tuned weights transfer approximately; per-cohort retune via Tier 2 typically improves accuracy by several millimetres. Similarly, the Tier 1.5 MLP is LEMON-trained; retraining per cohort is recommended for publication-grade accuracy on non-Polhemus data.

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

## Tier 1.5: MLP Residual Correction (opt-in)

A small learned residual correction can be layered on top of the geodesic prediction to reduce held-out error by ~70% on the LEMON cohort. The MLP learns systematic per-subject discrepancies between the geodesic's prediction and the true HTJ anatomy (notably the differential cap-standoff bias that a single rigid+scale Procrustes cannot absorb — see `Fiducial_Extrapolation_Exp/Results/c_anchor_augmentation_finding_20260419.md`).

**Architecture:** 1 hidden layer (15-input → 64 ReLU → 6-output), 1414 parameters, trained on 88 LEMON subjects' HTJ ground truth (2 QC-flagged subjects excluded). Inference is pure matmul (no ML-framework runtime dependency — numpy / MATLAB matrix ops only).

Enable by passing `mlp_correction=true` and a `cohort` tag:

```matlab
% MATLAB
[pL, pR, info] = geodesic_rescue(Cz, T7, T8, ...
    'cohort', 'LEMON_Polhemus_Adult', ...
    'mlp_correction', true);
% info.mlp_applied, info.mlp_correction_magnitude_L_mm/_R_mm surfaced
```

```python
# Python
from geodesic_rescue_py import apply_geodesic_rescue

raw_rescued = apply_geodesic_rescue(raw,
                                    cohort='LEMON_Polhemus_Adult',
                                    mlp_correction=True)
```

**Cohort-specific:** MLP weights are trained per-cohort. `LEMON_Polhemus_Adult` ships ready (4.76 mm held-out on LEMON). To deploy to a new cohort, retrain on pilot HTJ-tagged data using `Fiducial_Extrapolation_Exp/Scripts/c_Train_LEMON_MLP.py` as a template, then add an entry to `data/weight_zoo.json` pointing to the new weights file.

**Reproducibility:** the pre-MLP pure-geodesic state is preserved at git tag `v1.0.0-pure-geodesic`. Any paper citing the pre-V1.1 result should pin to that tag. Default path (without `mlp_correction`) is bit-identical to the tagged state — `c_Complete_Master_Scoreboard.m` wrapper parity test still passes 28/28 at 0.0000 mm.

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
