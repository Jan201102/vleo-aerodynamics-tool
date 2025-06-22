%% TORQUE CURVES
% Plot pitch torque over pitch angle for different models and geometries
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

%load environment data
run("environment_definitions.m");

%% Configuration
control_surface_angle__rad = 0.2; % 0 degrees for control surfaces

%% Load both geometries
% Load shuttlecock geometry from gmsh
bodies_shuttlecock = load_from_gmsh(energy_accommodation);
num_bodies_shuttlecock = 5;
rotation_face_index_shuttlecock = [2,3,4,5];

%% load lut data
lut_data = load_lut("aerodynamic_coefficients_panel_method.csv");

%% Pitch angle configuration
num_pitch_angles = 101;
pitch_angles__rad = linspace(-pi/6000, pi/6000, num_pitch_angles); % -30 to +30 degrees
axis_direction = [0; 1; 0]; % Rotation axis (y-axis for pitch)
torque_component = 2; % Component index (2 = y-component for pitch)

%% Calculate torque curves
% Initialize arrays for results
torque_shuttlecock = zeros(num_pitch_angles, 2); % 2 models

% Set control surface angles to 0 for shuttlecock geometry
bodies_rotation_angles_shuttlecock = zeros(1, num_bodies_shuttlecock);
bodies_rotation_angles_shuttlecock(rotation_face_index_shuttlecock) = control_surface_angle__rad;

% Calculate torque for each pitch angle and model (only shuttlecock)
for model = 1:2
    for i = 1:num_pitch_angles
        pitch_angle = pitch_angles__rad(i);
        
        % Calculate torque for shuttlecock geometry
        attitude_quaternion_BI = [cos(pitch_angle/2); sin(pitch_angle/2) * axis_direction];
        
        [~, torque_vec,~,~] = vleoAerodynamics(...
            attitude_quaternion_BI,...
            rotational_velocity_BI_B__rad_per_s,...
            velocity_I_I__m_per_s,...
            wind_velocity_I_I__m_per_s,...
            density__kg_per_m3,...
            temperature__K,...
            particles_mass__kg,...
            bodies_shuttlecock,...
            bodies_rotation_angles_shuttlecock,...
            temperature_ratio_method,...
            model,...
            lut_data);
        torque_shuttlecock(i, model) = torque_vec(torque_component);
    end
end

%% Calculate aerodynamic stiffness (tangent slopes) at pitch angle 0
derivation_method = 'central';
delta_alpha__rad = 1e-4;

% Calculate stiffness for shuttlecock geometry and models
stiffness_shuttlecock = zeros(2, 1);

for model = 1:2
    stiffness_shuttlecock(model) = calculate_aerodynamic_stiffness(derivation_method, ...
                                                                 delta_alpha__rad, ...
                                                                 axis_direction, ...
                                                                 torque_component, ...
                                                                 bodies_shuttlecock, ...
                                                                 bodies_rotation_angles_shuttlecock, ...
                                                                 model, ...
                                                                 lut_data);
end

% Calculate torque at pitch angle 0 for tangent line intercepts
zero_angle_idx = find(abs(pitch_angles__rad) == min(abs(pitch_angles__rad)), 1);
torque_at_zero_shuttlecock = torque_shuttlecock(zero_angle_idx, :);

% Calculate tangent lines
tangent_shuttlecock = zeros(num_pitch_angles, 2);

for model = 1:2
    tangent_shuttlecock(:, model) = torque_at_zero_shuttlecock(model) + stiffness_shuttlecock(model) * pitch_angles__rad;
end

%% Plot results
figure;
plot(rad2deg(pitch_angles__rad), torque_shuttlecock(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(rad2deg(pitch_angles__rad), torque_shuttlecock(:,2), 'r-', 'LineWidth', 2);

% Plot tangent lines
plot(rad2deg(pitch_angles__rad), tangent_shuttlecock(:,1), 'b:', 'LineWidth', 1.5);
plot(rad2deg(pitch_angles__rad), tangent_shuttlecock(:,2), 'r:', 'LineWidth', 1.5);

grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');
legend(...
    sprintf('Shuttlecock - Sentman \\alpha_E = %.2f', energy_accommodation), ...
    'Shuttlecock - IRS model', ...
    'Shuttlecock - Sentman Tangent', ...
    'Shuttlecock - IRS Tangent', ...
    'Location', 'best');
title(sprintf('Pitch Torque vs Pitch Angle with Tangents at 0° (Control Surface Angle = %.1f°) \\Delta_\\alpha = %.5f°', rad2deg(control_surface_angle__rad),rad2deg(delta_alpha__rad)));
hold off;
