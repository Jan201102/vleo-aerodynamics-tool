%% LOAD BOX TRUNCATED
%% This functions loads the box shaped satellite
%  from six separate .obj files for each side
function bodies = load_box_truncated(energy_accommodation_coefficient,surface_temp__K)
    arguments
        energy_accommodation_coefficient (1,1) double = 0.9; % Energy accommodation coefficient
        surface_temp__K (1,1) double = 300; % Surface temperature in Kelvin
    end

    import vleo_aerodynamics_core.*
    % Get absolute path of test folder
    [test_folder,~,~] = fileparts(mfilename("fullpath"));

    disp(test_folder)
    obj_files = string(fullfile(test_folder, '../obj', ...
                                 {'body_front.obj', ...
                                 'body_left.obj', ...
                                 'body_back.obj',...
                                 'body_right.obj',...
                                 'body_top.obj',...
                                 'body_bottom.obj'
                                 }));

    % Define rotation hinge points for each side (no rotation for box sides)
    rotation_hinge_points_CAD = zeros(3,6);

    % Define rotation directions for each side (no rotation for box sides)
    rotation_directions_CAD = ones(3,6);

    % Surface temperatures for all six sides
    surface_temperatures__K = num2cell(surface_temp__K*ones(1,6));

    % Surface energy accommodation coefficients for all six sides
    surface_energy_accommodation_coefficients = num2cell(energy_accommodation_coefficient*ones(1,6));

    % Use the same DCM and CoM as in load_box
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