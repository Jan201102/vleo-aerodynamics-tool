%% LOAD box
%% This functions loads the box shaped sattelite
%  from .obj files stored in the example folder
function bodies = load_box(energy_accommodation_coefficient,surface_temp__K)
    arguments
        energy_accommodation_coefficient (1,1) double = 0.9; % Energy accommodation coefficient
        surface_temp__K (1,1) double = 300; % Surface temperature in Kelvin
    end

    import vleo_aerodynamics_core.*
    % Get absolute path of test folder
    [test_folder,~,~] = fileparts(mfilename("fullpath"));

    disp(test_folder)
    obj_files = string(fullfile(test_folder, '../obj', 'body.obj'));
    rotation_hinge_points_CAD = zeros(3,1);

    rotation_directions_CAD = [1; 0; 0];

    surface_temperatures__K = num2cell(surface_temp__K);

    surface_energy_accommodation_coefficients = num2cell(energy_accommodation_coefficient);

    DCM_B_from_CAD = [0, -1, 0;...
                        -1, 0, 0; ...
                        0, 0, -1];

    CoM_CAD = [0; 2; 0];


    bodies = importMultipleBodies(obj_files, ...
        rotation_hinge_points_CAD, ...
        rotation_directions_CAD, ...
        surface_temperatures__K, ...
        surface_energy_accommodation_coefficients, ...
        DCM_B_from_CAD, ...
        CoM_CAD);
end