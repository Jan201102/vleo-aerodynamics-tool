%% LOAD SHUTTLECOCK WING
%% This function loads only the WingRight configuration from the shuttlecock satellite
function bodies = load_shuttlecock_wing()
    import vleo_aerodynamics_core.*
    
    % Get absolute path of test folder
    [test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    
    % Load only WingRight.obj (index 4 in the original array)
    object_files = string(fullfile(test_folder, 'obj', 'WingRight.obj'));
    
    % Use parameters for WingRight (index 4 from load_model)
    surface_temperatures__K = {300}; % 300*(3) for index 4 (0-based: 300*3)
    surface_energy_accommodation_coefficients = {0.9};
    
    % WingRight parameters (column 4 from load_model)
    rotation_hinge_points_CAD = [0.05; 0; 0];
    rotation_directions_CAD = [0; -1; 0];
    
    DCM_B_from_CAD = [0, 0, 1;...
                      0, 1, 0; ...
                      -1, 0, 0];
    CoM_CAD = [0; 0; 0.1];
    
    bodies = importMultipleBodies(object_files, ...
        rotation_hinge_points_CAD, ...
        rotation_directions_CAD, ...
        surface_temperatures__K, ...
        surface_energy_accommodation_coefficients, ...
        DCM_B_from_CAD, ...
        CoM_CAD);
end
