function bodies = load_unit_cube(energy_accommodation_coefficient, surface_temp__K)
%LOAD_UNIT_CUBE Load unit cube geometry for aerodynamic analysis
%
%   bodies = LOAD_UNIT_CUBE(energy_accommodation_coefficient, surface_temp__K)
%
%   Loads a unit cube geometry consisting of six separate faces from 
%   individual .obj files.
%
% Inputs:
%   energy_accommodation_coefficient - Double, energy accommodation coefficient (default: 0.9)
%   surface_temp__K                 - Double, surface temperature [K] (default: 300)
%
% Outputs:
%   bodies - Cell array containing body structures for each cube face

    arguments
        energy_accommodation_coefficient (1,1) double = 0.9
        surface_temp__K (1,1) double = 300
    end

    import vleo_aerodynamics_core.*
    
    %% Get File Paths
    [test_folder, ~, ~] = fileparts(mfilename('fullpath'));
    fprintf('Loading unit cube from: %s\n', test_folder);
    
    obj_files = string(fullfile(test_folder, '../obj/cube', ...
                                {'front.obj', ...
                                 'left.obj', ...
                                 'back.obj', ...
                                 'right.obj', ...
                                 'top.obj', ...
                                 'bottom.obj'}));

    %% Define Geometric Properties
    % Rotation hinge points for each face (no rotation for cube faces)
    rotation_hinge_points_CAD = zeros(3, 6);

    % Rotation directions for each face (no rotation for cube faces)
    rotation_directions_CAD = ones(3, 6);

    % Surface temperatures for all six faces
    surface_temperatures__K = num2cell(surface_temp__K*ones(1,6));

    % Surface energy accommodation coefficients for all six sides
    surface_energy_accommodation_coefficients = num2cell(energy_accommodation_coefficient*ones(1,6));

    % Use the same DCM and CoM as in load_box
    DCM_B_from_CAD = eye(3);

    CoM_CAD = [0; 0; 0];

    bodies = importMultipleBodies(obj_files, ...
        rotation_hinge_points_CAD, ...
        rotation_directions_CAD, ...
        surface_temperatures__K, ...
        surface_energy_accommodation_coefficients, ...
        DCM_B_from_CAD, ...
        CoM_CAD);
end