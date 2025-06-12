%% TORQUE CURVES
% Plot pitch torque over pitch angle for different models and geometries
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

%load environment data
run("environment_definitions.m");

%% Configuration
energy_accommodation = 0.9;
control_surface_angle = 0.0; % 0 degrees for control surfaces

%% Load both geometries
% Load shuttlecock geometry from gmsh
bodies_shuttlecock = load_from_gmsh(energy_accommodation);
num_bodies_shuttlecock = 5;
rotation_face_index_shuttlecock = [2,3,4,5];

% Load model geometry from obj files
bodies_model = load_model(energy_accommodation);
num_bodies_model = 5;
rotation_face_index_model = [2,3,4,5];

%% Load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
lut = fullfile(test_folder, 'aerodynamic_coefficients_panel_method.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% Pitch angle configuration
num_pitch_angles = 101;
pitch_angles__rad = linspace(-pi/6000, pi/6000, num_pitch_angles); % -30 to +30 degrees
axis_direction = [0; 1; 0]; % Rotation axis (y-axis for pitch)
torque_component = 2; % Component index (2 = y-component for pitch)

%% Calculate torque curves
% Initialize arrays for results
torque_shuttlecock = zeros(num_pitch_angles, 2); % 2 models
torque_model = zeros(num_pitch_angles, 2); % 2 models

% Set control surface angles to 0 for both geometries
bodies_rotation_angles_shuttlecock = zeros(1, num_bodies_shuttlecock);
bodies_rotation_angles_shuttlecock(rotation_face_index_shuttlecock) = control_surface_angle;

bodies_rotation_angles_model = zeros(1, num_bodies_model);
bodies_rotation_angles_model(rotation_face_index_model) = control_surface_angle;

% Calculate torque for each pitch angle and model
for model = 1:2
    for i = 1:num_pitch_angles
        pitch_angle = pitch_angles__rad(i);
        
        % Calculate torque for shuttlecock geometry
        attitude_quaternion_BI = [cos(pitch_angle/2); sin(pitch_angle/2) * axis_direction];
        
        [~, torque_vec] = vleoAerodynamics(...
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
            1,...
            lut);
        torque_shuttlecock(i, model) = torque_vec(torque_component);
        
        % Calculate torque for model geometry
        [~, torque_vec] = vleoAerodynamics(...
            attitude_quaternion_BI,...
            rotational_velocity_BI_B__rad_per_s,...
            velocity_I_I__m_per_s,...
            wind_velocity_I_I__m_per_s,...
            density__kg_per_m3,...
            temperature__K,...
            particles_mass__kg,...
            bodies_model,...
            bodies_rotation_angles_model,...
            temperature_ratio_method,...
            model,...
            1,...
            lut);
        torque_model(i, model) = torque_vec(torque_component);
    end
end

%% Calculate aerodynamic stiffness (tangent slopes) at pitch angle 0
derivation_method = 'central';
delta_alpha = 1e-4;

% Calculate stiffness for both geometries and models
stiffness_shuttlecock = zeros(2, 1);
stiffness_model = zeros(2, 1);

for model = 1:2
    stiffness_shuttlecock(model) = calculate_aerodynamic_stiffness(derivation_method, ...
                                                                 delta_alpha, ...
                                                                 axis_direction, ...
                                                                 torque_component, ...
                                                                 bodies_shuttlecock, ...
                                                                 bodies_rotation_angles_shuttlecock, ...
                                                                 model, ...
                                                                 lut);
    
    stiffness_model(model) = calculate_aerodynamic_stiffness(derivation_method, ...
                                                           delta_alpha, ...
                                                           axis_direction, ...
                                                           torque_component, ...
                                                           bodies_model, ...
                                                           bodies_rotation_angles_model, ...
                                                           model, ...
                                                           lut);
end

% Calculate torque at pitch angle 0 for tangent line intercepts
zero_angle_idx = find(abs(pitch_angles__rad) == min(abs(pitch_angles__rad)), 1);
torque_at_zero_shuttlecock = torque_shuttlecock(zero_angle_idx, :);
torque_at_zero_model = torque_model(zero_angle_idx, :);

% Calculate tangent lines
tangent_shuttlecock = zeros(num_pitch_angles, 2);
tangent_model = zeros(num_pitch_angles, 2);

for model = 1:2
    tangent_shuttlecock(:, model) = torque_at_zero_shuttlecock(model) + stiffness_shuttlecock(model) * pitch_angles__rad;
    tangent_model(:, model) = torque_at_zero_model(model) + stiffness_model(model) * pitch_angles__rad;
end

%% Plot results
figure;
plot(pitch_angles__rad, torque_shuttlecock(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(pitch_angles__rad, torque_shuttlecock(:,2), 'r-', 'LineWidth', 2);
plot(pitch_angles__rad, torque_model(:,1), 'b--', 'LineWidth', 2);
plot(pitch_angles__rad, torque_model(:,2), 'r--', 'LineWidth', 2);

% Plot tangent lines
plot(pitch_angles__rad, tangent_shuttlecock(:,1), 'b:', 'LineWidth', 1.5);
plot(pitch_angles__rad, tangent_shuttlecock(:,2), 'r:', 'LineWidth', 1.5);
plot(pitch_angles__rad, tangent_model(:,1), 'b-.', 'LineWidth', 1.5);
plot(pitch_angles__rad, tangent_model(:,2), 'r-.', 'LineWidth', 1.5);

grid on;
xlabel('Pitch Angle [rad]');
ylabel('Pitch Torque [Nm]');
legend(...
    sprintf('Shuttlecock - Sentman \\alpha_E = %.2f', energy_accommodation), ...
    'Shuttlecock - IRS model', ...
    sprintf('Load Model - Sentman \\alpha_E = %.2f', energy_accommodation), ...
    'Load Model - IRS model', ...
    'Shuttlecock - Sentman Tangent', ...
    'Shuttlecock - IRS Tangent', ...
    'Load Model - Sentman Tangent', ...
    'Load Model - IRS Tangent', ...
    'Location', 'best');
title(sprintf('Pitch Torque vs Pitch Angle with Tangents at 0° (Control Surface Angle = %.1f°)', control_surface_angle * 180/pi));

% Add text annotations showing stiffness values
text_x = max(pitch_angles__rad) * 0.7;
text_y_start = max([max(torque_shuttlecock(:)), max(torque_model(:))]) * 0.8;
text(text_x, text_y_start, sprintf('Stiffness Values [Nm/rad]:'), 'FontSize', 10, 'FontWeight', 'bold');
text(text_x, text_y_start*0.9, sprintf('Shuttlecock Sentman: %.3e', stiffness_shuttlecock(1)), 'FontSize', 9);
text(text_x, text_y_start*0.85, sprintf('Shuttlecock IRS: %.3e', stiffness_shuttlecock(2)), 'FontSize', 9);
text(text_x, text_y_start*0.8, sprintf('Load Model Sentman: %.3e', stiffness_model(1)), 'FontSize', 9);
text(text_x, text_y_start*0.75, sprintf('Load Model IRS: %.3e', stiffness_model(2)), 'FontSize', 9);

hold off;

%% Save figure
saveas(gcf, 'pitch_torque_curves_with_tangents.png');
saveas(gcf, 'pitch_torque_curves_with_tangents.eps', 'epsc');

%% Plot aerodynamic stiffness comparison
figure;
categories = {'Shuttlecock\nSentman', 'Shuttlecock\nIRS', 'Load Model\nSentman', 'Load Model\nIRS'};
stiffness_values = [stiffness_shuttlecock(1), stiffness_shuttlecock(2), stiffness_model(1), stiffness_model(2)];
colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980];

bar_plot = bar(1:4, stiffness_values, 'FaceColor', 'flat');
bar_plot.CData = colors;

% Add value labels on top of bars
for i = 1:4
    text(i, stiffness_values(i) + 0.05*max(abs(stiffness_values)), ...
         sprintf('%.3e', stiffness_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

grid on;
set(gca, 'XTickLabel', categories);
ylabel('Aerodynamic Stiffness [Nm/rad]');
title(sprintf('Aerodynamic Stiffness Comparison (Control Surface Angle = %.1f°)', control_surface_angle * 180/pi));

% Add legend
legend('Sentman Model', 'IRS Model', 'Location', 'best');

% Save stiffness comparison figure
saveas(gcf, 'aerodynamic_stiffness_comparison.png');
saveas(gcf, 'aerodynamic_stiffness_comparison.eps', 'epsc');

fprintf('Torque curves with tangents plotted and saved successfully!\n');
fprintf('Aerodynamic stiffness comparison plotted and saved successfully!\n');
fprintf('Aerodynamic stiffness values:\n');
fprintf('Shuttlecock - Sentman: %.6e Nm/rad\n', stiffness_shuttlecock(1));
fprintf('Shuttlecock - IRS: %.6e Nm/rad\n', stiffness_shuttlecock(2));
fprintf('Load Model - Sentman: %.6e Nm/rad\n', stiffness_model(1));
fprintf('Load Model - IRS: %.6e Nm/rad\n', stiffness_model(2));
