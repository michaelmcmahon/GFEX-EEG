# GFEX-EEG Reproduction Quickstart

This folder contains a fast end-to-end sanity check that the GFEX-EEG toolbox loads, runs, and produces deterministic output without requiring any external dataset download. It is intended primarily for **reviewers and first-time users** who want to verify the toolbox works on their machine in under a minute, before investing time in downloading the public EEG-MRI cohorts referenced in the paper.

## Files

| File | Purpose |
|---|---|
| `quickstart.py` | Python end-to-end demo + golden-value assertions |
| `quickstart.m` | MATLAB equivalent (same input, same expected output) |
| `expected_output.json` | Canonical predicted HTJ coordinates for the four demo modes |

## What the quickstart does

For each script (`quickstart.py` and `quickstart.m`), the demo:

1. Loads the canonical template Cz / T7 / T8 anchor positions (the same coordinates that ship in `matlab/core/template_anchors.m` and `python/geodesic_rescue_py/template_anchors.py`) and uses them as a synthetic subject input.
2. Runs the GFEX-EEG engine in four modes:
   - **Demo 1.** Tier 1 pure-geodesic prediction (no learned correction).
   - **Demo 2.** Tier 1.5 with the `LEMON_Polhemus_Adult` MLP correction.
   - **Demo 3.** Tier 1.5 with the `WH_Neuromag70_Adult` MLP correction.
   - **Demo 4.** Tier 1.5 with the `CapTrak_Adult` MLP correction.
3. Asserts each predicted LHJ / RHJ matches a hardcoded golden value within 0.01 mm Euclidean (the engine is deterministic; passing margin in practice is machine-precision, < 1e-10 m).

Expected runtime: **< 1 second**. Expected output on a clean install: **all four demos PASS**.

## Running it

**Python** (with `geodesic_rescue_py` either pip-installed or accessible via the toolbox source tree):

```bash
cd <toolbox-root>/reproduce_results
python quickstart.py
```

**MATLAB** (from the MATLAB prompt, with no special path setup required — the script self-registers):

```matlab
>> cd <toolbox-root>/reproduce_results
>> quickstart
```

Both scripts will print per-demo status, the predicted LHJ / RHJ coordinates, and the deviation from the golden value. A summary `PASS / FAIL` is printed at the end.

## What this verifies

A successful run confirms that:

- The MATLAB and Python core engines load without error.
- The ICBM152 scalp mesh (`data/ICBM152_scalp.mat`) is present and parses correctly.
- The cohort preset loader (`data/weight_zoo.json`) resolves all three production cohort tags (`LEMON_Polhemus_Adult`, `WH_Neuromag70_Adult`, `CapTrak_Adult`).
- All three production MLP weight files (`data/mlp/mlp_*.mat`) load without error.
- The full Procrustes alignment, Dijkstra geodesic walk, tangent extension, pivot-lock projection, and MLP residual correction pipeline runs end-to-end and produces deterministic, bit-identical output.

## What this does NOT verify

This script **does not** reproduce the per-cohort accuracy numbers reported in the GFEX-EEG manuscript (4.76 mm held-out HTJ error on LEMON, 8.16 mm on Wakeman-Henson, 8.52 mm on Shirley Ryan AbilityLab). Those numbers are accuracies against per-subject MRI ground truth, which requires the corresponding public EEG-MRI datasets:

- **LEMON**: distributed via [NITRC Fcon\_1000](http://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON.html).
- **Wakeman-Henson**: OpenNeuro `ds000117` (MRI) and `ds002718` (EEG); see the paper's Methods for the biological-fingerprint re-pairing procedure.
- **Shirley Ryan AbilityLab TMS-EEG-MRI-fMRI-DWI**: OpenNeuro `ds004024`.

Subject-specific HTJ tags performed in Brainstorm are available from the corresponding author on reasonable request (see the manuscript's Data Availability statement).

## How to extend the quickstart

If you have your own per-subject Cz / T7 / T8 anchor coordinates (in any unit; the engine auto-detects and normalises), you can substitute them for the template anchors at the top of either script. The golden-value assertions will then no longer apply (those are pinned to the template input), but the rest of the pipeline will run identically and report predicted LHJ / RHJ in the same RAS-metre output frame.

For full retuning on a new hardware class (e.g., a non-Polhemus / non-Neuromag / non-CapTrak cohort), see the Tier 2 metaheuristic tuning section of the toolbox README and the `matlab/tuning/geodesic_fno_tuner.m` / `geodesic_rescue_py.GeodesicTuner` entry points.

## Reporting issues

If `quickstart.py` or `quickstart.m` reports a `FAIL` on a clean install, please file an issue at <https://github.com/michaelmcmahon/GFEX-EEG/issues> with:

- The full stdout output of the failing script.
- Your operating system + Python version (or MATLAB version).
- The toolbox release tag you have checked out (`git describe --tags`).
