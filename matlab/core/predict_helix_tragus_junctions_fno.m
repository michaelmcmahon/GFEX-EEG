function [LHJ_mni, RHJ_mni, iCz, iT7, iT8] = predict_helix_tragus_junctions_fno(Cz_t, T7_t, T8_t, mesh_path, rho, beta, D_L_sub, D_R_sub, D_standard)
%PREDICT_HELIX_TRAGUS_JUNCTIONS_FNO  V2.2 geodesic walk on the ICBM152 scalp mesh.
%
%   [LHJ_mni, RHJ_mni, iCz, iT7, iT8] = PREDICT_HELIX_TRAGUS_JUNCTIONS_FNO( ...
%       Cz_t, T7_t, T8_t, mesh_path, rho, beta, D_L_sub, D_R_sub, D_standard)
%
%   Performs the V2.2 ("projection paradox fix") two-arm geodesic descent
%   from the Cz seed across the scalp mesh, returns the predicted left
%   and right HTJ in template (MNI) space, along with the closest-vertex
%   indices used for the Procrustes mapping back to subject space.
%
%   This file is LOCKED — do not modify without an explicit re-tune of
%   the rho / beta parameters and full cross-platform parity re-test.
%
% ==============================================================================
%   GFEX-EEG TOOLBOX — Geodesic walk (MATLAB, V2.2, LOCKED)
% ------------------------------------------------------------------------------
%   Authors:
%     Michael McMahon  (ORCID: 0000-0002-5266-3194)
%     Michael Schukat  (ORCID: 0000-0002-6908-6100)
%     Enda Barrett     (ORCID: 0000-0002-9876-8717)
%     University of Galway, Galway, Ireland
%
%   Repository : https://github.com/michaelmcmahon/GFEX-EEG
%   Issues     : https://github.com/michaelmcmahon/GFEX-EEG/issues
%
%   CITATION (please cite both)
%     [Software] McMahon, M., Schukat, M., & Barrett, E. (2026).
%                GFEX-EEG Toolbox [Software].
%                Zenodo. https://doi.org/10.5281/zenodo.20580899
%     [Paper]    McMahon, M., Schukat, M., & Barrett, E. (Submitted).
%                GFEX-EEG: Geodesic recovery of anatomical fiducials for
%                MRI-free EEG source imaging.
%     See also CITATION.cff in the repository root (machine-readable).
%
%   VERSION
%     For the current package version, call gfex_version().
%
%   LICENSE
%     SPDX-License-Identifier: MIT
%     SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway
%     See LICENSE in the repository root.
% ==============================================================================

    if nargin < 9, D_standard = 0.1388; end
    mesh = load(mesh_path);
    V = mesh.Vertices; F = mesh.Faces;
    % 1. Build Mesh Graph
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = unique(sort(E, 2), 'rows');
    dist = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2, 2));
    G = graph(E(:,1), E(:,2), dist);

    % 2. Manifold Registration (Capture Indices for Procrustes)
    [~, iCz] = min(sum(bsxfun(@minus, V, Cz_t(:)').^2, 2));
    [~, iT7] = min(sum(bsxfun(@minus, V, T7_t(:)').^2, 2));
    [~, iT8] = min(sum(bsxfun(@minus, V, T8_t(:)').^2, 2));
    
    % 3. Individualized Ratio using Subject Scalars
    Ratio_L = rho * ((D_L_sub / D_standard) ^ beta);
    Ratio_R = rho * ((D_R_sub / D_standard) ^ beta);

    % 4. Geodesic Path Traversal (MNI Space)
    path_L = shortestpath(G, iCz, iT7);
    path_R = shortestpath(G, iCz, iT8);
    
    d_path_L = sum(sqrt(sum(diff(V(path_L,:)).^2, 2)));
    d_path_R = sum(sqrt(sum(diff(V(path_R,:)).^2, 2)));
    
    target_s_L = d_path_L * (1 + Ratio_L);
    target_s_R = d_path_R * (1 + Ratio_R);

    % 5. Free-Space Tangent Continuation
    LHJ_mni = trace_geodesic_distance(V, path_L, target_s_L, d_path_L);
    RHJ_mni = trace_geodesic_distance(V, path_R, target_s_R, d_path_R);
end

function final_pt = trace_geodesic_distance(vertices, path, target_s, path_s)
    extra_dist = target_s - path_s;
    idx_start = max(1, length(path) - 5);
    tangent_vec = vertices(path(end),:) - vertices(path(idx_start),:);
    tangent_vec = tangent_vec / norm(tangent_vec);
    final_pt = vertices(path(end),:) + (tangent_vec * extra_dist);
end
