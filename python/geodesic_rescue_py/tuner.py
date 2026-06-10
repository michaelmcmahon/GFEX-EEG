"""GFEX-EEG: Far/Near metaheuristic hyperparameter tuner (Python)

Re-tunes the geodesic descent ratio (rho) and non-linear scaling exponent
(beta) for non-standard EEG hardware classes. Two-phase optimization:

    Far : global search via differential_evolution over
          rho in [0.001, 0.1], beta in [0.5, 2.5]
    Near: local refinement via Nelder-Mead simplex with TolX ~1e-4

Loss is the radial-telescope metric (V37.0) which decouples direction
from magnitude and is robust to outliers in the held-out pool.

Usage
-----
    from geodesic_rescue_py import GeodesicTuner
    tuner = GeodesicTuner(training_data=[...])
    rho_opt, beta_opt = tuner.tune()

Notes
-----
Bit-identical loss-function output with the MATLAB tuner
(geodesic_fno_tuner.m).

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
Current version: see ``geodesic_rescue_py.__version__`` (single source
of truth in ``__init__.py``).

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

import numpy as np
from scipy.optimize import minimize, differential_evolution
from .core import GeodesicRescue

class GeodesicTuner:
    """
    FNO Metaheuristic Optimizer for Geodesic Rescue parameters.
    Used to derive custom rho and beta for specialized hardware or cohorts.
    """
    def __init__(self, training_data, mesh_path=None):
        """
        training_data: list of dicts with keys:
            'Cz', 'T7', 'T8': [3] np.array (Meters)
            'gtLPA', 'gtRPA': [3] np.array (Meters ground truth)
        """
        self.data = training_data
        self.rescuer = GeodesicRescue(mesh_path=mesh_path)

    def _objective_function(self, theta):
        rho, beta = theta
        if rho < 0 or beta < 0:
            return 1e6
            
        errs = []
        for s in self.data:
            pL, pR = self.rescuer.rescue(
                s['Cz'], s['T7'], s['T8'], 
                rho=rho, beta=beta
            )
            
            # V37.0 Radial Telescope Logic
            tm_pred = (s['T7'] + s['T8']) / 2.0
            tm_gt = (s['gtLPA'] + s['gtRPA']) / 2.0
            
            V_pred_L = pL - tm_pred
            V_gt_L = s['gtLPA'] - tm_gt
            V_pred_R = pR - tm_pred
            V_gt_R = s['gtRPA'] - tm_gt
            
            norm_V_pred_L = np.linalg.norm(V_pred_L)
            norm_V_pred_R = np.linalg.norm(V_pred_R)
            
            if norm_V_pred_L > 0:
                V_telescope_L = (V_pred_L / norm_V_pred_L) * np.linalg.norm(V_gt_L)
            else:
                V_telescope_L = V_pred_L
                
            if norm_V_pred_R > 0:
                V_telescope_R = (V_pred_R / norm_V_pred_R) * np.linalg.norm(V_gt_R)
            else:
                V_telescope_R = V_pred_R
                
            err_L = np.linalg.norm(V_telescope_L - V_gt_L) * 1000
            err_R = np.linalg.norm(V_telescope_R - V_gt_R) * 1000
            
            # Calculate mean Euclidean error in mm
            err = (err_L + err_R) / 2.0
            errs.append(err)
            
        return np.mean(errs)

    def tune(self, max_iter_far=20, max_iter_near=50):
        """
        Executes dual-phase FNO optimization.
        """
        print("\n>>> INITIATING FNO PHASE 1: FAR SEARCH (Global Exploration) <<<")
        
        # Bounds for rho and beta
        bounds = [(0.001, 0.100), (0.5, 2.5)]
        
        # Differential Evolution for robust global search (Phase 1: FAR)
        res_far = differential_evolution(
            self._objective_function, 
            bounds, 
            maxiter=max_iter_far,
            popsize=10,
            disp=True
        )
        
        theta_far = res_far.x
        print(f"   [FAR] Best Initial Basin: Rho={theta_far[0]:.4f}, Beta={theta_far[1]:.4f} (Residual: {res_far.fun:.2f} mm)")

        print("\n>>> INITIATING FNO PHASE 2: NEAR SEARCH (Local Refinement) <<<")
        
        # L-BFGS-B or Nelder-Mead for local refinement (Phase 2: NEAR)
        res_near = minimize(
            self._objective_function,
            theta_far,
            method='Nelder-Mead',
            options={'maxiter': max_iter_near, 'disp': True}
        )
        
        theta_opt = res_near.x
        final_residual = res_near.fun
        
        print("\n--- FNO TUNING COMPLETE ---")
        print(f"Optimized Rho:  {theta_opt[0]:.6f}")
        print(f"Optimized Beta: {theta_opt[1]:.6f}")
        print(f"Final Mean Residual: {final_residual:.4f} mm\n")
        
        return {
            'optimal_rho': theta_opt[0],
            'optimal_beta': theta_opt[1],
            'final_residual_mm': final_residual,
            'success': res_near.success
        }
