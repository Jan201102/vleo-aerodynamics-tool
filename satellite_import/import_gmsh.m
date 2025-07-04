% ========================================================================
% GMSH_IMPORT: Import a .m-file generated by GMSH (containing a 'msh' struct)
%              and convert it into a structured format for body-wise processing.
%
% Usage:
%   satellite = gmsh_import('C:/path/to/mesh/full_satellite.m');
%
% The input file must define a variable 'msh' with fields:
%   - msh.POS:     Node positions (n_nodes × 3)
%   - msh.TRIANGLES: Element definitions (n_elements × 4), where
%                    columns 1:3 = node indices of each triangle
%                    column 4   = body ID
%
% Output:
%   body_data struct with:
%     - vertices_B:  cell array of size {n_bodies×1}, each cell is 3×3×n_triangles
%     - centroids_B: cell array {n_bodies×1}, each 3×n matrix of triangle centroids
%     - normals_B:   cell array {n_bodies×1}, each 3×n matrix of normal vectors
%     - areas_B:     cell array {n_bodies×1}, each 1×n vector of triangle areas
% ========================================================================

function body_data = import_gmsh(relative_path)

    % Save current working directory to restore it later
    old_dir = pwd;

    % Separate file path, name and extension
    [filepath, filename, ext] = fileparts(relative_path);

    % If a folder is specified, switch into it
    if ~isempty(filepath)
        cd(filepath);
    end

    clear msh;  % Ensure no leftover msh variable from workspace

    % Run the GMSH-exported .m file (should define 'msh')
    run(fullfile(filepath, [filename, ext]));

    % Check if the file exists
    if ~exist(fullfile(filepath, [filename, ext]), 'file')
        error('File not found: %s', fullfile(filepath, [filename, ext]));
    end

    % Validate msh struct existence
    if ~exist('msh', 'var')
        error('The GMSH-exported file does not define a "msh" structure.');
    end

    % Get all unique body IDs from the 4th column of the triangle data
    body_ids = unique(msh.TRIANGLES(:,4));

    % Extract vertex coordinates for each body
    % For each body:
    % - collect all triangle node indices (3 per triangle)
    % - use these to extract XYZ coordinates from msh.POS
    % - reshape to a 3×3×n array (3 nodes × 3 coordinates × n triangles)
    vertices_B = arrayfun(@(b) ...
        reshape(msh.POS(reshape(msh.TRIANGLES(msh.TRIANGLES(:,4) == b, 1:3)', [], 1), :)', ...
                [3, 3, nnz(msh.TRIANGLES(:,4) == b)]), ...
        body_ids, 'UniformOutput', false);

    % Compute triangle centroids per body (mean of each triangle's 3 vertices)
    centroids_B = cellfun(@(V) ...
        reshape(mean(V, 2), 3, []), ...
        vertices_B, ...
        'UniformOutput', false);

    % Compute normal vectors for each triangle: cross(v2−v1, v3−v1)
    normals_B = cellfun(@(V) ...
        reshape(cross(V(:,2,:) - V(:,1,:), V(:,3,:) - V(:,1,:)), 3, []), ...
        vertices_B, ...
        'UniformOutput', false);

    % Compute triangle areas: 0.5 * norm of each normal vector
    areas_B = cellfun(@(N) ...
        0.5 * sqrt(sum(N.^2, 1)), ...
        normals_B, ...
        'UniformOutput', false);

    % Return to original working directory
    cd(old_dir);

    % Create output structure
    body_data = struct('vertices_B', vertices_B, ...
                       'centroids_B', centroids_B, ...
                       'normals_B', normals_B, ...
                       'areas_B', areas_B);

end

