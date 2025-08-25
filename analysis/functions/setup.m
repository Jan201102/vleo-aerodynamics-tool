function [bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut] = setup(lut_file, geometry_type, energy_accommodation)
    arguments
        lut_file string
        geometry_type string
        energy_accommodation double = -1
    end

    import vleo_aerodynamics_core.*;

    % Constants
    environment_definitions = struct();
    %DSMC parameter:
    % T_i=934K
    % T_w = 300K
    % u = 7800 m/s
    % n = 4.698e14 1/m^3
    % partikel_mass = 16U (Atomarer auserstoff) -> rho = 1.2482e-11 kg/m^3
    environment_definitions.rotational_velocity_BI_B__rad_per_s = 0;
    environment_definitions.velocity_I_I__m_per_s = 7800 * [1;0;0];
    environment_definitions.wind_velocity_I_I__m_per_s = zeros(3,1);
    environment_definitions.n = 4.698e14;
    environment_definitions.temperature__K = 934;
    environment_definitions.surface_temperature__K = 300;
    environment_definitions.particles_mass__kg = 16 * 1.6605390689252e-27;
    environment_definitions.density__kg_per_m3 = environment_definitions.particles_mass__kg * environment_definitions.n;
    environment_definitions.temperature_ratio_method = 1;
    if energy_accommodation < 0
        environment_definitions.energy_accommodation = (7.5e-17 * environment_definitions.n * environment_definitions.temperature__K)/(1+7.5e-17 * environment_definitions.n * environment_definitions.temperature__K);
    else
        environment_definitions.energy_accommodation = energy_accommodation;
    end

    %load lut based on file ending
    if endsWith(lut_file, '.mat')
        lut = load(lut_file);
    else
        [parent_folder,~,~] = fileparts(mfilename('fullpath'));
        parent_folder = fullfile(parent_folder, '..');
        display(parent_folder)
        lut_path = fullfile(parent_folder,lut_file);
        if ~isfile(lut_path)
            error("Look-up table file not found. Please check the path: %s", lut_path);
        end
        lut = readmatrix(lut_path);
        lut = griddedInterpolant(lut(:,1), lut(:,2:3),"spline");
    end

    %load geometry
    [bodies, num_bodies, rotation_face_index, x_label] = load_geometry(geometry_type,...
                                                                       environment_definitions.energy_accommodation,...
                                                                       environment_definitions.surface_temperature__K);

end