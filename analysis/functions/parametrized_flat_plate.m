function bodies = parametrized_flat_plate(x_dim, y_dim, cog_x, cog_y, two_sided,energy_accommodation_coefficient,surface_temp__K)
    % CREATE_FLAT_PLATE Creates a flat plate body structure
    %
    % Inputs:
    %   x_dim - dimension in x direction (width)
    %   y_dim - dimension in y direction (height)
    %   cog_x - x coordinate of center of gravity in CAD system
    %   cog_y - y coordinate of center of gravity in CAD system
    %   two_sided - logical, true for two-sided plate, false for one-sided (optional, default: true)
    %
    % Output:
    %   bodies - 1x1 cell containing body structure in body frame (centered at COG)
    arguments
        x_dim (1,1) double {mustBePositive} = 1.0; % Width of the plate
        y_dim (1,1) double {mustBePositive} = 1.0; % Height of the plate
        cog_x (1,1) double = 0.0; % Center of gravity x coordinate in CAD system
        cog_y (1,1) double = 0.0; % Center of gravity y coordinate in CAD system
        two_sided (1,1) logical = true; % Two-sided plate flag
        energy_accommodation_coefficient (1,1) double = 0.9; % Energy accommodation coefficient
        surface_temp__K (1,1) double = 300; % Surface temperature in Kelvin
    end
    % Default to two-sided if not specified
    if nargin < 5
        two_sided = true;
    end
    
    % Half dimensions for easier vertex calculation
    half_x = x_dim / 2;
    half_y = y_dim / 2;
    
    % Create body structure
    body = struct();
    
    % Determine number of triangles based on plate type
    num_triangles = two_sided * 4 + ~two_sided * 2;
    
    % Create vertices in CAD system (plate symmetric about origin)
    % Then transform to body frame by subtracting COG coordinates
    body.vertices_B = zeros(3,3,num_triangles);
    
    % Bottom face triangles (normal pointing in -z direction)
    % CAD vertices transformed to body frame
    body.vertices_B(:,:,1) = [-half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y, -half_y - cog_y,  half_y - cog_y;
                              -1e-10, -1e-10, -1e-10];
    
    body.vertices_B(:,:,2) = [ half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y,  half_y - cog_y,  half_y - cog_y;
                              -1e-10, -1e-10, -1e-10];
    
    if two_sided
        % Top face triangles (normal pointing in +z direction)
        body.vertices_B(:,:,3) = [-half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                                  -half_y - cog_y, -half_y - cog_y,  half_y - cog_y;
                                   1e-10,  1e-10,  1e-10];
        
        body.vertices_B(:,:,4) = [ half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                                  -half_y - cog_y,  half_y - cog_y,  half_y - cog_y;
                                   1e-10,  1e-10,  1e-10];
    end
    
    % Calculate centroids
    body.centroids_B = squeeze(mean(body.vertices_B, 2));
    
    % Define normals (pointing outward)
    if two_sided
        body.normals_B = [0 0 0 0; 0 0 0 0; -1 -1 1 1];
    else
        body.normals_B = [0 0; 0 0; -1 -1];
    end
    
    % Calculate areas (each triangle is half the total plate area)
    triangle_area = (x_dim * y_dim) / 2;
    body.areas = repmat(triangle_area, 1, num_triangles);
    
    % Rotation hinge point at body frame origin (COG in CAD system)
    body.rotation_hinge_point_B = [0; 0; 0];
    
    % Rotation direction (y-axis for pitch)
    body.rotation_direction_B = [0; -1; 0];
    
    % Default material properties
    body.temperatures__K = repmat(surface_temp__K, num_triangles, 1);
    body.energy_accommodation_coefficients = repmat(energy_accommodation_coefficient, num_triangles, 1);
    
    % Return as 1x1 cell
    bodies = cell(1,1);
    bodies{1} = body;
end
