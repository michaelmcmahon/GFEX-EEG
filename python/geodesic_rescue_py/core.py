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
import scipy.io
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
import os

def procrustes_analysis(P_subj, P_temp):
    """
    Standard Procrustes analysis with scaling (Umeyama algorithm).
    Maps P_temp to P_subj.
    """
    mu_subj = P_subj.mean(0)
    mu_temp = P_temp.mean(0)
    
    P_subj_c = P_subj - mu_subj
    P_temp_c = P_temp - mu_temp
    
    sig_subj = np.sum(P_subj_c**2) / P_subj.shape[0]
    sig_temp = np.sum(P_temp_c**2) / P_temp.shape[0]
    
    # SVD
    C = (P_subj_c.T @ P_temp_c) / P_subj.shape[0]
    U, S, Vt = np.linalg.svd(C)
    
    # Rotation
    R = U @ Vt
    if np.linalg.det(R) < 0:
        # Reflection fix
        S_fix = np.eye(U.shape[1])
        S_fix[-1, -1] = -1
        R = U @ S_fix @ Vt
        
    # Scale
    b = np.sum(S) / sig_temp if sig_temp > 0 else 1.0
    
    # Translation
    c = mu_subj - b * (mu_temp @ R.T)
    
    return b, R, c

def predict_helix_tragus_junctions_fno(Cz_t, T7_t, T8_t, vertices, faces, rho, beta, D_L_sub, D_R_sub, D_standard=0.140):
    """
    Python translation of the Geodesic Extrapolation Algorithm.
    """
    # 1. Build Mesh Graph
    num_verts = vertices.shape[0]
    
    # Edges from faces
    edges = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    
    # Distances
    dist = np.linalg.norm(vertices[edges[:, 0]] - vertices[edges[:, 1]], axis=1)
    
    # Sparse adjacency matrix for Dijkstra
    adj = csr_matrix((dist, (edges[:, 0], edges[:, 1])), shape=(num_verts, num_verts))
    adj = adj + adj.T  # Symmetric
    
    # 2. Find closest indices
    def find_closest(pt):
        return np.argmin(np.sum((vertices - pt)**2, axis=1))
    
    iCz = find_closest(Cz_t)
    iT7 = find_closest(T7_t)
    iT8 = find_closest(T8_t)
    
    # 3. Individualized Ratio
    Ratio_L = rho * ((D_L_sub / D_standard) ** beta)
    Ratio_R = rho * ((D_R_sub / D_standard) ** beta)
    
    # 4. Geodesic Path Traversal
    # Scipy returns the distances and predecessors
    dist_L, pred_L = shortest_path(adj, directed=False, indices=iCz, return_predecessors=True)
    dist_R, pred_R = shortest_path(adj, directed=False, indices=iCz, return_predecessors=True)
    
    def get_path(predecessors, start, end):
        path = [end]
        while path[-1] != start:
            path.append(predecessors[path[-1]])
            if path[-1] == -9999: # No path
                return []
        return path[::-1]
    
    path_L = get_path(pred_L, iCz, iT7)
    path_R = get_path(pred_R, iCz, iT8)
    
    def path_length(path_indices):
        if not path_indices: return 0
        return np.sum(np.linalg.norm(np.diff(vertices[path_indices], axis=0), axis=1))
    
    d_path_L = path_length(path_L)
    d_path_R = path_length(path_R)
    
    target_s_L = d_path_L * (1 + Ratio_L)
    target_s_R = d_path_R * (1 + Ratio_R)
    
    def trace_geodesic_distance(path, target_s, path_s):
        extra_dist = target_s - path_s
        idx_start = max(0, len(path) - 6)
        tangent_vec = vertices[path[-1]] - vertices[path[idx_start]]
        norm_val = np.linalg.norm(tangent_vec)
        if norm_val > 0:
            tangent_vec = tangent_vec / norm_val
        return vertices[path[-1]] + (tangent_vec * extra_dist)

    LHJ_mni = trace_geodesic_distance(path_L, target_s_L, d_path_L)
    RHJ_mni = trace_geodesic_distance(path_R, target_s_R, d_path_R)
    
    return LHJ_mni, RHJ_mni, iCz, iT7, iT8

class GeodesicRescue:
    def __init__(self, mesh_path=None):
        if mesh_path is None:
            # Default to bundled ICBM152
            base_dir = os.path.dirname(os.path.abspath(__file__))
            mesh_path = os.path.join(base_dir, 'data', 'ICBM152_scalp.mat')
        
        data = scipy.io.loadmat(mesh_path)
        self.vertices = data['Vertices']
        self.faces = data['Faces'] - 1 # 0-indexed
        
    def rescue(self, Cz_m, T7_m, T8_m, rho=0.0054, beta=1.1885):
        """
        Main rescue function. Expects coordinates in meters.
        Returns predicted LHJ, RHJ in meters.
        """
        # Template Anchors (Golden Kite in Meters)
        T_Cz = np.array([0.011240, 0.025921, 0.141134])
        T_T7 = np.array([-0.089174, -0.001327, -0.006348])
        T_T8 = np.array([0.096880, -0.014286, -0.005819])
        
        D_L_sub = np.linalg.norm(Cz_m - T7_m)
        D_R_sub = np.linalg.norm(Cz_m - T8_m)
        
        # 1. Execute Core Engine (MNI Space Prediction)
        Lp_mni, Rp_mni, iCz, iT7, iT8 = predict_helix_tragus_junctions_fno(
            T_Cz, T_T7, T_T8, self.vertices, self.faces, rho, beta, D_L_sub, D_R_sub
        )
        
        # 2. Temporal Pivot Mapping (Subject Space)
        P_temp = np.vstack([self.vertices[iCz], self.vertices[iT7], self.vertices[iT8]])
        P_subj = np.vstack([Cz_m, T7_m, T8_m])
        
        b, R, c = procrustes_analysis(P_subj, P_temp)
        
        # Mapping using pivot-lock logic from MATLAB
        t_pivot = (T_T7 + T_T8) / 2
        s_pivot = (T7_m + T8_m) / 2
        
        # MATLAB: pLHJ = tr.b * (Lp_mni(1,:) - t_pivot) * tr.T + s_pivot;
        # Python: pLHJ = b * (Lp_mni - t_pivot) @ R.T + s_pivot
        pLHJ = b * (Lp_mni - t_pivot) @ R.T + s_pivot
        pRHJ = b * (Rp_mni - t_pivot) @ R.T + s_pivot
        
        return pLHJ, pRHJ
