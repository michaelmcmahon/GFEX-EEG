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
    MATLAB-exact procrustes: matches procrustes(P_subj, P_temp, 'scaling', true, 'reflection', false).
    Returns (b, T, c) such that P_subj ~= b * P_temp @ T + c. Apply as b * Y @ T + c (NOT @ T.T).
    """
    X = P_subj
    Y = P_temp

    muX = X.mean(0)
    muY = Y.mean(0)
    X0 = X - muX
    Y0 = Y - muY

    normX = np.sqrt(np.sum(X0**2))
    normY = np.sqrt(np.sum(Y0**2))

    X0n = X0 / normX
    Y0n = Y0 / normY

    A = X0n.T @ Y0n
    U, S, Vt = np.linalg.svd(A)
    L = U
    M = Vt.T
    D = S.copy()

    T = M @ L.T
    if np.linalg.det(T) < 0:
        M = M.copy()
        M[:, -1] = -M[:, -1]
        D[-1] = -D[-1]
        T = M @ L.T

    traceTA = np.sum(D)
    b = traceTA * normX / normY
    c = muX - b * (muY @ T)

    return b, T, c

def predict_helix_tragus_junctions_fno(Cz_t, T7_t, T8_t, vertices, faces, rho, beta, D_L_sub, D_R_sub, D_standard=0.1388):
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
    # Hard defaults (mirror weight_zoo.json 'default' alias)
    _DEFAULT_RHO        = 0.248383
    _DEFAULT_BETA       = 0.235926
    _DEFAULT_D_STANDARD = 0.1388

    def __init__(self, mesh_path=None, rho=None, beta=None, parity='RAS',
                 D_standard=None, cohort=None):
        # Resolve weights with precedence: explicit kwarg > cohort preset > hard default
        if cohort is not None:
            from .weight_zoo import load_cohort_preset
            preset = load_cohort_preset(cohort)
            if rho is None:        rho = preset['rho']
            if beta is None:       beta = preset['beta']
            if D_standard is None: D_standard = preset['D_standard']
        if rho is None:        rho        = self._DEFAULT_RHO
        if beta is None:       beta       = self._DEFAULT_BETA
        if D_standard is None: D_standard = self._DEFAULT_D_STANDARD

        is_default_mesh = False
        if mesh_path is None:
            base_dir = os.path.dirname(os.path.abspath(__file__))
            mesh_path = os.path.join(base_dir, 'data', 'ICBM152_scalp.mat')
            is_default_mesh = True

        data = scipy.io.loadmat(mesh_path)
        self.vertices = data['Vertices']
        # ICBM152_scalp.mat is already in RAS (T_T8 X=+0.097 = right = RAS).
        # Do NOT apply any axis conversion here.

        self.faces = data['Faces'] - 1 # 0-indexed
        self.rho = rho
        self.beta = beta
        self.parity = parity
        self.D_standard = D_standard
        self.cohort = cohort

    def rescue(self, Cz_m, T7_m, T8_m, rho=None, beta=None, D_standard=None,
               cohort=None):
        """
        Main rescue function. Expects coordinates in meters.
        Returns predicted LHJ, RHJ in meters.
        """
        # Per-call cohort override (rare; usually set at __init__ time)
        if cohort is not None:
            from .weight_zoo import load_cohort_preset
            preset = load_cohort_preset(cohort)
            if rho is None:         rho = preset['rho']
            if beta is None:        beta = preset['beta']
            if D_standard is None:  D_standard = preset['D_standard']

        rho_val = self.rho if rho is None else rho
        beta_val = self.beta if beta is None else beta
        D_std_val = self.D_standard if D_standard is None else D_standard
        
        # 1. PARITY RESTORE (Input ALS -> RAS)
        Cz_proc = Cz_m.copy()
        T7_proc = T7_m.copy()
        T8_proc = T8_m.copy()
        
        if self.parity.upper() in ['ALS', 'BIDS-BRAINSTORM']:
            Cz_proc = np.array([-Cz_m[1], Cz_m[0], Cz_m[2]])
            T7_proc = np.array([-T7_m[1], T7_m[0], T7_m[2]])
            T8_proc = np.array([-T8_m[1], T8_m[0], T8_m[2]])

        # 2. TEMPLATE ANCHORS (Exact MATLAB V22.2 Parity)
        T_Cz = np.array([0.011240, 0.025921, 0.141134])
        T_T7 = np.array([-0.089174, -0.001327, -0.006348])
        T_T8 = np.array([0.096880, -0.014286, -0.005819])

        # 3. INDIVIDUALIZED SCALING (Hardware Standoff)
        D_L_sub = np.linalg.norm(Cz_proc - T7_proc)
        D_R_sub = np.linalg.norm(Cz_proc - T8_proc)

        # 4. CORE ENGINE (Geodesic Trace in MNI)
        Lp_mni, Rp_mni, iCz, iT7, iT8 = predict_helix_tragus_junctions_fno(
            T_Cz, T_T7, T_T8, self.vertices, self.faces, rho_val, beta_val, D_L_sub, D_R_sub, D_std_val
        )
        
        # 5. PIVOT-LOCK PROJECTION (MNI -> Subject RAS)
        # Matches MATLAB geodesic_rescue.m: tr.b*(Lp_mni-t_pivot)*tr.T + s_pivot
        P_temp = np.vstack([self.vertices[iCz], self.vertices[iT7], self.vertices[iT8]])
        P_subj = np.vstack([Cz_proc, T7_proc, T8_proc])

        b, T, _ = procrustes_analysis(P_subj, P_temp)

        t_pivot = (T_T7 + T_T8) / 2
        s_pivot = (T7_proc + T8_proc) / 2
        pLHJ = b * (Lp_mni - t_pivot) @ T + s_pivot
        pRHJ = b * (Rp_mni - t_pivot) @ T + s_pivot
        
        return pLHJ, pRHJ
