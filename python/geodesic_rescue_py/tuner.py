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
        if rho <= 0 or beta <= 0:
            return 1e6
            
        errs = []
        for s in self.data:
            pL, pR = self.rescuer.rescue(
                s['Cz'], s['T7'], s['T8'], 
                rho=rho, beta=beta
            )
            # Calculate mean Euclidean error in mm
            err = (np.linalg.norm(pL - s['gtLPA']) + np.linalg.norm(pR - s['gtRPA'])) / 2 * 1000
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
