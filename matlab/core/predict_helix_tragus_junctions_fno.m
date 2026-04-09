% /*******************************************************************************
% * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
% * Version: 1.0.0
% * Repository: https://github.com/michaelmcmahon/GFEX-EEG
% * License:  MIT License
% * Authors: Michael McMahon / University of Galway
% * DOI: [If available]
% * Date: 2026
% *
% * [License Text or link to License file]
% *******************************************************************************/

% =========================================================================
% GEODESIC EXTRAPOLATION ALGORITHM (V2.2 - PROJECTION PARADOX FIX)
% Objective: Predict HTJ and return mesh indices for rigorous mapping.
% =========================================================================
function [LHJ_mni, RHJ_mni, iCz, iT7, iT8] = predict_helix_tragus_junctions_fno(Cz_t, T7_t, T8_t, mesh_path, rho, beta, D_L_sub, D_R_sub, D_standard)

    if nargin < 9, D_standard = 0.140; end
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
