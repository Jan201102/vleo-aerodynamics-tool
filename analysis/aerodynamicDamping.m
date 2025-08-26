%% AERODYNAMIC DAMPING ANALYSIS
% Calculate the aerodynamic damping for a satellite for different control 
% surface configurations using both Sentman and IRS models
import vleo_aerodynamics_core.*
clear;

%% Setup Environment, Geometry, and LUT
geometry_type = 'shuttlecock';  % Options: 'plate' or 'shuttlecock'
lut_file = 'aerodynamic_coefficients_panel_method_poly.mat';
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);

%% Derivation Configuration
derivation_method = 'central'; 
delta_q_rad_per_s = 1e2;          % Step size for numerical differentiation [rad/s]
axis_direction = [0; 1; 0];       % Rotation axis (y-axis for pitch)
torque_component = 2;             % Component index (2 = y-component for pitch)

%% Calculate Aerodynamic Damping
num_angles = 100;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aero_damping = zeros(num_angles, 2);

for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        aero_damping(i, model) = calculate_aerodynamic_damping(...
            derivation_method, ...
            delta_q_rad_per_s, ...
            axis_direction, ...
            torque_component, ...
            bodies, ...
            bodies_rotation_angles__rad, ...
            model, ...
            lut_data, ...
            environment_definitions);
        
        fprintf('Calculated point %d of %d\n', i + (model-1)*num_angles, 2*num_angles);
    end
end

% Save results
save('aerodynamic_damping_shuttlecock.mat', 'aero_damping', 'control_surface_angles__rad');
%% Plot Results
figure;
plot(rad2deg(control_surface_angles__rad), aero_damping(:,1), 'b');
hold on;
plot(rad2deg(control_surface_angles__rad), aero_damping(:,2), 'r');
grid on;
xlabel(x_label);
ylabel('Aerodynamic Damping [Nm/(rad/s)]');
legend(sprintf('Sentman \\alpha_E = %.4f', environment_definitions.energy_accommodation), ...
       'IRS model', 'Location', 'best');
title('Aerodynamic Damping vs Control Surface Angle');
hold off;

%% Save figure as PNG and EPS
matlab2tikz(sprintf('aerodynamic_damping_%s.tex', geometry_type));
