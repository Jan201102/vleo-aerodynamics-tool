import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

% Import constants from environment:definitions.m
run('environment_definitions.m');

%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'cl_cd_cVAE_A01_flat_and_bird.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

sattelite = load_model();

showBodies(sattelite, [0,0/4,pi/4,pi/4,pi/4], 0.75, 0.25);

%%
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles, 2);
aerodynamic_torque_B_B__Nm = aerodynamic_force_B__N;
attitude_quaternion_BI = [1; 0; 0; 0];
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1,5);
        bodies_rotation_angles__rad(4) = current_angle;
        [aerodynamic_force_B__N(:,i,model), ...
            aerodynamic_torque_B_B__Nm(:,i,model)] = ...
            vleoAerodynamics(attitude_quaternion_BI, ...
                             rotational_velocity_BI_B__rad_per_s, ...
                             velocity_I_I__m_per_s, ...
                             wind_velocity_I_I__m_per_s, ...
                             density__kg_per_m3, ...
                             temperature__K, ... 
                             particles_mass__kg, ...
                             sattelite, ...                                                       
                             bodies_rotation_angles__rad, ...
                             temperature_ratio_method,...
                             model,...
                             lut);
    end
end

figure;
title('Aerodynamic Torque efficiency');
hold on;
grid on;
plot(aerodynamic_torque_B_B__Nm(2,:,1), -aerodynamic_force_B__N(1,:,1), 'b', 'DisplayName', 'Sentman');
plot(aerodynamic_torque_B_B__Nm(2,:,2), -aerodynamic_force_B__N(1,:,2), 'r', 'DisplayName', 'IRS model');
xlabel('pitch Torque [Nm]');
ylabel('Drag [N]');
legend();