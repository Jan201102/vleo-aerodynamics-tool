function [bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut] = ...
         setup(lut_file, geometry_type, energy_accommodation)
    %SETUP Initialize environment, geometry, and lookup table for aerodynamic analysis
    %
    % Inputs:
    %   lut_file              - String, path to lookup table file
    %   geometry_type         - String, geometry type ('plate', 'shuttlecock', etc.)
    %   energy_accommodation  - Double, optional energy accommodation coefficient
    %
    % Outputs:
    %   bodies                - Body structure array
    %   num_bodies           - Number of bodies
    %   rotation_face_index  - Index of rotatable face
    %   x_label              - Label for x-axis in plots
    %   environment_definitions - Environment parameter structure
    %   lut                  - Lookup table data
    
    arguments
        lut_file string
        geometry_type string
        energy_accommodation double = -1
    end

    import vleo_aerodynamics_core.*;

    %% Environment Definitions
    environment_definitions = struct();
    
    % DSMC parameters based on VLEO conditions:
    % T_i = 934K (gas temperature)
    % T_w = 300K (wall/surface temperature)  
    % u = 7800 m/s (orbital velocity)
    % n = 4.698e14 1/m³ (number density)
    % Particle mass = 16U (atomic oxygen)
    % Resulting density: rho = 1.2482e-11 kg/m³
    
    environment_definitions.rotational_velocity_BI_B__rad_per_s = 0;
    environment_definitions.velocity_I_I__m_per_s = 7800 * [1; 0; 0];
    environment_definitions.wind_velocity_I_I__m_per_s = zeros(3, 1);
    environment_definitions.n = 4.698e14;
    environment_definitions.temperature__K = 934;
    environment_definitions.surface_temperature__K = 300;
    environment_definitions.particles_mass__kg = 16 * 1.6605390689252e-27;
    environment_definitions.density__kg_per_m3 = ...
        environment_definitions.particles_mass__kg * environment_definitions.n;
    environment_definitions.temperature_ratio_method = 1;
    
    % Calculate energy accommodation coefficient
    if energy_accommodation < 0
        environment_definitions.energy_accommodation = ...
            (7.5e-17 * environment_definitions.n * environment_definitions.temperature__K) / ...
            (1 + 7.5e-17 * environment_definitions.n * environment_definitions.temperature__K);
    else
        environment_definitions.energy_accommodation = energy_accommodation;
    end
    %% Load Lookup Table
    % Load LUT based on file extension
    if endsWith(lut_file, '.mat')
        lut = load(lut_file);
    else
        [parent_folder, ~, ~] = fileparts(mfilename('fullpath'));
        parent_folder = fullfile(parent_folder, '..');
        fprintf('Parent folder: %s\n', parent_folder);
        lut_path = fullfile(parent_folder, lut_file);
        
        if ~isfile(lut_path)
            error("Look-up table file not found. Please check the path: %s", lut_path);
        end
        
        lut = readmatrix(lut_path);
        lut = griddedInterpolant(lut(:,1), lut(:,2:3), "spline");
    end

    %% Load Geometry
    [bodies, num_bodies, rotation_face_index, x_label] = load_geometry(...
        geometry_type, ...
        environment_definitions.energy_accommodation, ...
        environment_definitions.surface_temperature__K);

end