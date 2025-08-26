function bodies = parametrized_flat_plate(x_dim, y_dim, cog, two_sided, energy_accommodation_coefficient, surface_temp__K)
%PARAMETRIZED_FLAT_PLATE Create a parametrized flat plate body structure
%
%   bodies = PARAMETRIZED_FLAT_PLATE(x_dim, y_dim, cog, two_sided, energy_accommodation_coefficient, surface_temp__K)
%
%   Creates a flat plate body structure with specified dimensions and properties.
%   The plate can be configured as one-sided or two-sided for aerodynamic analysis.
%
% Inputs:
%   x_dim                           - Double, plate width [m]
%   y_dim                           - Double, plate height [m] 
%   cog                             - 3x1 double, center of gravity coordinates [m]
%   two_sided                       - Logical, true for two-sided plate (default: true)
%   energy_accommodation_coefficient - Double, energy accommodation coefficient (default: 0.9)
%   surface_temp__K                 - Double, surface temperature [K] (default: 300)
%
% Outputs:
%   bodies - 1x1 cell array containing body structure in body frame (centered at COG)

    arguments
        x_dim (1,1) double {mustBePositive} = 1.0
        y_dim (1,1) double {mustBePositive} = 1.0
        cog (3,1) double = [0; 0; 0]
        two_sided (1,1) logical = true
        energy_accommodation_coefficient (1,1) double = 0.9
        surface_temp__K (1,1) double = 300
    end
    
    % Half dimensions for easier vertex calculation
    half_x = x_dim / 2;
    half_y = y_dim / 2;
    
    %% Create Body Structure
    body = struct();
    
    % Determine number of triangles based on plate type
    num_triangles = two_sided * 4 + ~two_sided * 2;
    
    % Create vertices in CAD system (plate symmetric about origin)
    % Then transform to body frame by subtracting COG coordinates
    body.vertices_B = zeros(3,3,num_triangles);
    
    % Bottom face triangles (normal pointing in -z direction)
    % CAD vertices transformed to body frame
    body.vertices_B(:,:,1) = [-half_x,  half_x, -half_x;
                              -half_y, -half_y,  half_y;
                              -1e-10, -1e-10, -1e-10] - cog ;
    
    body.vertices_B(:,:,2) = [ half_x ,  half_x , -half_x ;
                              -half_y ,  half_y ,  half_y ;
                              -1e-10, -1e-10, -1e-10]- cog ;
    
    if two_sided
        % Top face triangles (normal pointing in +z direction)
        body.vertices_B(:,:,3) = [-half_x ,  half_x , -half_x ;
                                  -half_y , -half_y ,  half_y ;
                                   1e-10,  1e-10,  1e-10]- cog ;
        
        body.vertices_B(:,:,4) = [ half_x ,  half_x , -half_x ;
                                  -half_y ,  half_y ,  half_y ;
                                   1e-10,  1e-10,  1e-10]- cog;
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
