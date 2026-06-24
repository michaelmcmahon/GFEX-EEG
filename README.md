# GFEX-EEG Toolbox (V1.1.6)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20580898.svg)](https://doi.org/10.5281/zenodo.20580898)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Hardware-Agnostic Spatial Extrapolation Engine for EEG Fiducials**

GFEX-EEG solves the "MRI-Free, Fiducial-Free" paradox for high-density EEG systems (like EGI and CapTrak) and recovers floating origins from legacy datasets (like Wakeman-Henson). It uses a Far & Near Optimization (FNO) Metaheuristic paired with a Procrustes-Dijkstra Manifold Engine to extrapolate **Helix-Tragus Junction (HTJ)** coordinates — a geometrically unambiguous anatomical landmark — using only the positions of the Cz, T7, and T8 electrodes. For downstream compatibility, predicted HTJ coordinates are written into the standard `LPA`/`RPA` slots of MNE-Python Raw / EEGLAB / FieldTrip / Brainstorm data structures; existing source-imaging pipelines consume them as anatomical fiducials without modification.

**V1.1.6 refines the README wording and fixes a `setup.py` packaging bug** so `pip install .` succeeds from a clean clone. **V1.1.5 adds the MIT LICENSE file** at the repo root and the CITATION.cff metadata needed for the Zenodo deposit and GitHub's "Cite this repository" button. **V1.1.4 added a reviewer-friendly end-to-end sanity check** at [`reproduce_results/`](./reproduce_results) — `quickstart.py` and `quickstart.m` exercise all four production prediction modes (pure geodesic + three cohort MLPs) using shipped template anchors and golden-value assertions; no external dataset download required, runs in under a second. **V1.1.3 fixed the CapTrak_Adult MLP packaging and weights format** so all three production cohort MLPs (`LEMON_Polhemus_Adult` 4.76 mm, `WH_Neuromag70_Adult` 8.16 mm, `CapTrak_Adult` 8.52 mm held-out) are runnable through both the MATLAB and Python inference paths. **V1.1.0 introduced the optional Tier 1.5 residual-correction MLP**; opt in via the `mlp_correction=true` kwarg. The default pure-geodesic path is preserved bit-identically and remains tagged at `v1.0.0-pure-geodesic` for reproducibility of pre-V1.1 results. See the Tier 1.5 section below.

## Quickstart

To verify the toolbox loads and runs correctly on your machine, without downloading any external dataset:

```bash
cd reproduce_results
python quickstart.py
```

```matlab
>> cd reproduce_results
>> quickstart
```

Both scripts run all four production prediction modes on shipped template anchors and assert against hardcoded golden values. Expected output on a clean install: **all 4 demos PASS** in well under a second. See [`reproduce_results/README.md`](./reproduce_results/README.md) for details and instructions on extending the quickstart to your own per-subject anchor data.

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
├── reproduce_results/             <- Quickstart sanity check (Python + MATLAB)
└── README.md
```

## Forensic Data Warning: The "Z=0 Coordinate Flatline" and "Raw Voxel Trap"
When extracting ground truth from legacy neuroinformatics pipelines, **always ensure your data hasn't been structurally flattened or grid-locked.** 

During our analysis of the LEMON dataset, we found that certain coordinate-system export conventions mathematically coerce the 3D Polhemus coordinates onto the `Z=0` plane, effectively **erasing the subject's physiological head tilt**. Default MRI coordinates are also frequently stored as **unscaled, un-tilted Voxel Indices**. 

**The Fix:** Always extract True World Space coordinates (`cs_convert` to 'world' in Brainstorm) to restore the subject's physiological tilt and 1.0mm scaling before running any FNO Tuner derivations. Failure to do so will result in an Orthogonal Parity Collision.

## Tier 1: "Black-Box" Wrappers (For Standard Adult Cohorts)

For 95% of use cases, researchers can rely on the LEMON-tuned weights (`rho = 0.248383`, `beta = 0.235926`, tuned against manually-annotated HTJ ground truth on N=10 LEMON subjects, held-out mean error 16.82 mm on N=10 independent subjects). The wrappers automatically apply the **Hyper-Scale Catch** to neutralize coordinate-system double-scaling traps and execute **RAS-to-ALS** axis transposition.

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

**Currently shipped production entries (held-out HTJ error against MRI ground truth, Tier 1.5 MLP-corrected):**

- `LEMON_Polhemus_Adult` — default; Polhemus per-subject digitization; **4.76 mm** (Leipzig LEMON cohort, `ds000221`, N=10 held-out)
- `WH_Neuromag70_Adult` — Neuromag 70-ch cap; **8.16 mm** (Wakeman-Henson cohort, `ds000117`/`ds002718`, N=5 held-out)
- `CapTrak_Adult` — BrainProducts CapTrak digitization; **8.52 mm** (Shirley Ryan AbilityLab TMS-EEG-MRI-fMRI-DWI cohort, `ds004024`, N=13 leave-one-out CV)

**Pending community contribution:** `EGI_HydroCel_256_Adult`, `Neuroscan_SynAmps_Adult` — entries are reserved in the zoo with `status: pending_tune`. Researchers with MRI-bearing cohorts on those setups can run the Tier 2 FNO tuner and submit a PR adding a production entry. See `data/weight_zoo.json` `contribution` block for the submission checklist.

### Validation & Accuracy

Held-out HTJ error against per-subject manually-tagged MRI ground truth, across three independent EEG-MRI cohorts spanning Polhemus, Neuromag, and BrainProducts CapTrak digitization hardware:

| Cohort | Hardware | N held-out | Tier 1 (pure geodesic) | Tier 1.5 (+ MLP) | Improvement |
|---|---|---|---|---|---|
| LEMON (`ds000221`) | Polhemus | 10 | 16.82 mm (SD 4.03) | **4.76 mm** (SD 2.06) | 71.7% |
| Wakeman-Henson (`ds000117`/`ds002718`) | Neuromag 70-ch | 5 | 21.60 mm (SD 6.44) | **8.16 mm** (SD 1.55) | 62.2% |
| Shirley Ryan AbilityLab (`ds004024`) | BrainProducts CapTrak | 13 (LOO-CV) | 18.36 mm (SD 6.31) | **8.52 mm** (SD 5.71) | 53.6% |

The LEMON 4.76 mm Tier 1.5 result is statistically indistinguishable from a measured ~5 mm inter-modality annotation floor (Polhemus stylus vs MRI tagging) across both LEMON and Wakeman-Henson cohorts, indicating that further algorithmic refinement on this manifold is precluded by human annotation noise rather than by method limitation. The Tier 1.5 result also survives a 4-test data-leak diligence battery (subject separation; label-shuffle permutation test with 3× degradation on scrambled labels; alternative-holdout rotation at 5.29 ± 0.62 mm across 5 random 10-subject splits; FNO-seen vs unseen within training at 0.15 mm gap). Full report: `Fiducial_Extrapolation_Exp/Results/c_MLP_Diligence_20260420.json`.

**Deploying to hardware outside the three shipped cohort presets.** The production presets — `LEMON_Polhemus_Adult` (Polhemus), `WH_Neuromag70_Adult` (Neuromag 70-ch), and `CapTrak_Adult` (BrainProducts CapTrak) — cover the dominant adult digitization classes. For other hardware (EGI HydroCel, Neuroscan, dense pediatric caps, etc.), the pure-geodesic Tier 1 path with the LEMON-tuned defaults transfers approximately and is often within several millimetres of post-retune accuracy. For publication-grade results on non-shipped hardware classes, run a Tier 2 FNO retune of ρ/β on a small pilot HTJ-tagged subset and (ideally) a Tier 1.5 MLP retrain on the same; submit the resulting preset back to the zoo via the `data/weight_zoo.json` `contribution` block.

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

**Cohort-specific:** MLP weights are trained per-cohort. Three cohort MLPs ship ready: `LEMON_Polhemus_Adult` (4.76 mm held-out, Polhemus), `WH_Neuromag70_Adult` (8.16 mm, Neuromag 70-ch), and `CapTrak_Adult` (8.52 mm, BrainProducts CapTrak). To deploy to a new cohort, retrain on pilot HTJ-tagged data using `Fiducial_Extrapolation_Exp/Scripts/c_Train_LEMON_MLP.py` as a template, then add an entry to `data/weight_zoo.json` pointing to the new weights file.

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

## Citation

If you use this software, please cite both the toolbox (via the Zenodo DOI) and the accompanying manuscript when it is published. The CITATION.cff file at the repository root provides ready-to-export APA / BibTeX / EndNote formats via GitHub's "Cite this repository" button.

**Toolbox archive:** [10.5281/zenodo.20580898](https://doi.org/10.5281/zenodo.20580898) (concept DOI — always resolves to the latest version; current release v1.1.6, 2026-06-24)

**Manuscript:** McMahon, M., Schukat, M., & Barrett, E. "GFEX-EEG: Geodesic recovery of anatomical fiducials for MRI-free EEG source imaging." (Submitted for publication.)

---

**Authors:** Michael McMahon, Michael Schukat, Enda Barrett — University of Galway  
**Repository:** [https://github.com/michaelmcmahon/GFEX-EEG](https://github.com/michaelmcmahon/GFEX-EEG)  
**License:** MIT License  
**Released:** 2026-06-24 (v1.1.6)

