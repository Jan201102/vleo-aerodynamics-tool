function bodies = parametrized_flat_plate(x_dim, y_dim, cog_x, cog_y)
    % CREATE_FLAT_PLATE Creates a flat plate body structure
    %
    % Inputs:
    %   x_dim - dimension in x direction (width)
    %   y_dim - dimension in y direction (height)
    %   cog_x - x coordinate of center of gravity in CAD system
    %   cog_y - y coordinate of center of gravity in CAD system
    %
    % Output:
    %   bodies - 1x1 cell containing body structure in body frame (centered at COG)
    
    % Half dimensions for easier vertex calculation
    half_x = x_dim / 2;
    half_y = y_dim / 2;
    
    % Create body structure
    body = struct();
    
    % Create vertices in CAD system (plate symmetric about origin)
    % Then transform to body frame by subtracting COG coordinates
    body.vertices_B = zeros(3,3,4);
    
    % Bottom face triangles (normal pointing in -z direction)
    % CAD vertices transformed to body frame
    body.vertices_B(:,:,1) = [-half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y, -half_y - cog_y,  half_y - cog_y;
                              -1e-10, -1e-10, -1e-10];
    
    body.vertices_B(:,:,2) = [ half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y,  half_y - cog_y,  half_y - cog_y;
                              -1e-10, -1e-10, -1e-10];
    
    % Top face triangles (normal pointing in +z direction)
    body.vertices_B(:,:,3) = [-half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y, -half_y - cog_y,  half_y - cog_y;
                               1e-10,  1e-10,  1e-10];
    
    body.vertices_B(:,:,4) = [ half_x - cog_x,  half_x - cog_x, -half_x - cog_x;
                              -half_y - cog_y,  half_y - cog_y,  half_y - cog_y;
                               1e-10,  1e-10,  1e-10];
    
    % Calculate centroids
    body.centroids_B = squeeze(mean(body.vertices_B, 2));
    
    % Define normals (pointing outward)
    body.normals_B = [0 0 0 0; 0 0 0 0; -1 -1 1 1];
    
    % Calculate areas (each triangle is half the total plate area)
    triangle_area = (x_dim * y_dim) / 2;
    body.areas = [triangle_area, triangle_area, triangle_area, triangle_area];
    
    % Rotation hinge point at body frame origin (COG in CAD system)
    body.rotation_hinge_point_B = [0; 0; 0];
    
    % Rotation direction (y-axis for pitch)
    body.rotation_direction_B = [0; -1; 0];
    
    % Default material properties
    body.temperatures__K = [300; 300; 300; 300];
    body.energy_accommodation_coefficients = [0.9; 0.9; 0.9; 0.9];
    
    % Return as 1x1 cell
    bodies = cell(1,1);
    bodies{1} = body;
end
