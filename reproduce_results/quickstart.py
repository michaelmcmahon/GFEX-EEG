"""
GFEX-EEG quickstart (Python)
============================
Sanity-check the toolbox end-to-end without downloading any external dataset.

Uses the shipped template Cz/T7/T8 anchors (canonical RAS-metre positions on
the ICBM152 scalp manifold) as a synthetic subject input, then exercises:

  1. Pure geodesic prediction        (Tier 1)
  2. LEMON_Polhemus_Adult MLP        (Tier 1.5)
  3. WH_Neuromag70_Adult MLP         (Tier 1.5)
  4. CapTrak_Adult MLP               (Tier 1.5)

For each mode the predicted Helix-Tragus Junction coordinates are compared
against hardcoded golden values (machine-precision tolerance: 1e-5 m = 0.01
millimetre). All four assertions should PASS on a clean install.

This script demonstrates that the toolbox loads, runs, and produces
deterministic bit-identical output across runs. It does NOT reproduce the
per-cohort accuracy numbers (4.76 / 8.16 / 8.52 mm) reported in the paper:
those additionally require downloading the public LEMON / Wakeman-Henson /
ds004024 datasets per the manuscript's Data Availability statement.

Usage:
    cd <toolbox-root>/reproduce_results
    python quickstart.py

Expected runtime: <1 second.
"""

import json
from pathlib import Path
import numpy as np

# Make the toolbox importable when running from the reproduce_results folder
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from geodesic_rescue_py import GeodesicRescue


# ---------------------------------------------------------------------------
# Synthetic subject input: the shipped template Cz/T7/T8 anchors
# ---------------------------------------------------------------------------
Cz = np.array([0.011240, 0.025921, 0.141134])
T7 = np.array([-0.089174, -0.001327, -0.006348])
T8 = np.array([0.096880, -0.014286, -0.005819])

TOL = 1e-5  # 0.01 mm tolerance for golden-value comparison

# Load expected golden values
expected_path = Path(__file__).parent / "expected_output.json"
with open(expected_path) as f:
    expected = json.load(f)


def check(label, lhj, rhj, expected_key):
    golden = expected[expected_key]
    exp_lhj = np.array(golden["LHJ"])
    exp_rhj = np.array(golden["RHJ"])
    delta_l = np.linalg.norm(lhj - exp_lhj)
    delta_r = np.linalg.norm(rhj - exp_rhj)
    pass_l = delta_l < TOL
    pass_r = delta_r < TOL
    status = "PASS" if (pass_l and pass_r) else "FAIL"
    print(f"  {label:30s}  delta L = {delta_l*1000:.6f} mm, delta R = {delta_r*1000:.6f} mm   [{status}]")
    print(f"    LHJ predicted: [{lhj[0]:+.6f}, {lhj[1]:+.6f}, {lhj[2]:+.6f}] m")
    print(f"    RHJ predicted: [{rhj[0]:+.6f}, {rhj[1]:+.6f}, {rhj[2]:+.6f}] m")
    assert pass_l and pass_r, f"{label} failed golden-value check"


print("GFEX-EEG quickstart (Python) -- end-to-end sanity check")
print("=" * 72)
print()
print("Input anchors (RAS metres):")
print(f"  Cz: [{Cz[0]:+.6f}, {Cz[1]:+.6f}, {Cz[2]:+.6f}]")
print(f"  T7: [{T7[0]:+.6f}, {T7[1]:+.6f}, {T7[2]:+.6f}]")
print(f"  T8: [{T8[0]:+.6f}, {T8[1]:+.6f}, {T8[2]:+.6f}]")
print()

gr = GeodesicRescue(cohort="LEMON_Polhemus_Adult")

print("Demo 1 -- Pure geodesic (Tier 1)")
lhj, rhj = gr.rescue(Cz, T7, T8)
check("LEMON cohort, pure geodesic", lhj, rhj, "pure_geodesic")
print()

print("Demo 2 -- LEMON_Polhemus_Adult MLP correction (Tier 1.5)")
lhj, rhj = gr.rescue(Cz, T7, T8, cohort="LEMON_Polhemus_Adult", mlp_correction=True)
check("LEMON_Polhemus_Adult + MLP", lhj, rhj, "lemon_mlp")
print()

print("Demo 3 -- WH_Neuromag70_Adult MLP correction (Tier 1.5)")
lhj, rhj = gr.rescue(Cz, T7, T8, cohort="WH_Neuromag70_Adult", mlp_correction=True)
check("WH_Neuromag70_Adult + MLP", lhj, rhj, "wh_mlp")
print()

print("Demo 4 -- CapTrak_Adult MLP correction (Tier 1.5)")
lhj, rhj = gr.rescue(Cz, T7, T8, cohort="CapTrak_Adult", mlp_correction=True)
check("CapTrak_Adult + MLP", lhj, rhj, "captrak_mlp")
print()

print("=" * 72)
print("All 4 demos PASSED. Toolbox is functioning end-to-end.")
print()
print("Note: this verifies deterministic engine behaviour against golden")
print("values. Reproducing the per-cohort accuracy numbers (4.76 mm LEMON,")
print("8.16 mm Wakeman-Henson, 8.52 mm CapTrak) reported in the manuscript")
print("additionally requires downloading the public datasets named in the")
print("paper's Data Availability statement.")
