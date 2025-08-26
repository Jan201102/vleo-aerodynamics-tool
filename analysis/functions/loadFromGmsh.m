function bodies = load_from_gmsh(energy_accommodation_coefficient, surface_temp__K)
%LOAD_FROM_GMSH Load shuttlecock geometry from GMSH-generated files
%
%   bodies = LOAD_FROM_GMSH(energy_accommodation_coefficient, surface_temp__K)
%
%   Loads the shuttlecock satellite configuration with main body and
%   four movable control surfaces from GMSH-generated .obj files.
%
% Inputs:
%   energy_accommodation_coefficient - Double, energy accommodation coefficient (default: 0.9)
%   surface_temp__K                 - Double, surface temperature [K] (default: 300)
%
% Outputs:
%   bodies - Cell array containing body structures for main body and control surfaces

    arguments
        energy_accommodation_coefficient (1,1) double = 0.9
        surface_temp__K (1,1) double = 300
    end

    import vleo_aerodynamics_core.*
    
    %% Get File Path
    [test_folder, ~, ~] = fileparts(mfilename('fullpath'));
    gmsh_file = fullfile(test_folder, '../obj', 'Shuttlecock.m');

    %% Surface Properties
    surface_temperatures__K = num2cell(surface_temp__K * ones(1, 5));
    surface_energy_accommodation_coefficients = ...
        num2cell(energy_accommodation_coefficient * ones(1, 5));

    %% Rotation Configuration
    % Rotation hinge points for each body [m]
    rotation_hinge_points_CAD = [0, -0.15, -0.15, -0.15, -0.15;
                                 0,  0.05,  0.00, -0.05,  0.00;
                                 0,  0.00, -0.05,  0.00,  0.05];

    % Rotation axis directions for each body
    rotation_directions_CAD = [0,  0,  0,  0,  0;
                               1,  0, -1,  0,  1;
                               0, -1,  0,  1,  0];

    %% Coordinate System Transformation
    DCM_B_from_CAD = eye(3);
    CoM_CAD = [0; 0; 0.0];

    %% Load Bodies
    bodies = importMultipleBodies(gmsh_file, ...
        rotation_hinge_points_CAD, ...
        rotation_directions_CAD, ...
        surface_temperatures__K, ...
        surface_energy_accommodation_coefficients, ...
        DCM_B_from_CAD, ...
        CoM_CAD);
end