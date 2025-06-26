%% LOAD MODEL
%% This functions loads the shuttle cock configuration sattelite
%  from .obj files stored in the example folder
function bodies = load_from_gmsh(energy_accommodation_coefficient,surface_temp__K)
    arguments
        energy_accommodation_coefficient (1,1) double = 0.9; % Energy accommodation coefficient
        surface_temp__K (1,1) double = 300; % Surface temperature in Kelvin
    end

    import vleo_aerodynamics_core.*
    % Get absolute path of test folder
    %[test_folder,~,~] = fileparts(mfilename("fullpath"));
    [test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    gmsh_file = fullfile(test_folder, 'obj', 'Shuttlecock_copy_BA.m');

    surface_temperatures__K = num2cell(surface_temp__K*ones(1,5));
    surface_energy_accommodation_coefficients = num2cell(energy_accommodation_coefficient*ones(1,5));

    rotation_hinge_points_CAD = [0,-0.15,-0.15,-0.15,-0.15;...
                                0,0.05,0,-0.05,0;...
                                0,0,-0.05,0,0.05];

    rotation_directions_CAD = [0,0,0,0,0;...
                                1,0,-1,0,1;...
                                0,-1,0,1,0];

    DCM_B_from_CAD = eye(3);
    CoM_CAD = [0; 0; 0.0];

    bodies = importMultipleBodies(gmsh_file,...
    rotation_hinge_points_CAD, ...
    rotation_directions_CAD, ...
    surface_temperatures__K, ...
    surface_energy_accommodation_coefficients, ...
    DCM_B_from_CAD, ...
    CoM_CAD);
end