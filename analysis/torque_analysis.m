%% TORQUE ANALYSIS
% Analyze aerodynamic torques for unit cube geometry
% Compare Sentman and IRS models across different pitch angles
import vleo_aerodynamics_core.*;

%% Setup Environment, Geometry, and LUT
geometry_type = 'unit_cube';
lut_file = "aerodynamic_coefficients_panel_method_poly.mat";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = ...
    setup(lut_file, geometry_type);
bodies_label = ["bottom", "left", "top", "right", "back", "front"];

%% Calculate Torque Curves
% Set control surface angles to 0 for unit cube geometry
bodies_rotation_angles = zeros(num_bodies);
num_pitch_angles = 201;
pitch_angles__rad = linspace(0, 2*pi, num_pitch_angles);
torques = zeros(num_pitch_angles, 2, 3, num_bodies);  % [angles, models, xyz, bodies]
total_torques = zeros(num_pitch_angles, 2, 3);        % [angles, models, xyz]
axis_direction = [0; 1; 0];  % Rotation axis (y-axis for pitch)

% Calculate torque for each pitch angle and model
for model = 1:2
    for i = 1:num_pitch_angles
        pitch_angle = pitch_angles__rad(i);
        
        % Calculate attitude quaternion for current pitch angle
        attitude_quaternion_BI = [cos(pitch_angle/2); sin(pitch_angle/2) * axis_direction];
        
        [~, total_torque_vec, ~, torque_vec] = vleoAerodynamics(...
            attitude_quaternion_BI, ...
            environment_definitions.rotational_velocity_BI_B__rad_per_s, ...
            environment_definitions.velocity_I_I__m_per_s, ...
            environment_definitions.wind_velocity_I_I__m_per_s, ...
            environment_definitions.density__kg_per_m3, ...
            environment_definitions.temperature__K, ...
            environment_definitions.particles_mass__kg, ...
            bodies, ...
            bodies_rotation_angles, ...
            environment_definitions.temperature_ratio_method, ...
            model, ...
            lut_data);
        
        torques(i, model, :, :) = torque_vec;
        total_torques(i, model, :) = total_torque_vec';
    end
end
%% Plot Pitch Torques for Each Body
figure;
for i = 1:num_bodies
    subplot(3, 2, i);
    hold on;
    
    % Plot torques for each model (y-component = pitch)
    plot(rad2deg(pitch_angles__rad), torques(:,1,2,i), 'b-', 'LineWidth', 2);
    plot(rad2deg(pitch_angles__rad), torques(:,2,2,i), 'r-', 'LineWidth', 2);
    
    grid on;
    xlabel('Pitch Angle [°]');
    ylabel('Pitch Torque [Nm]');
    title(sprintf('%s Face', bodies_label(i)));
    
    if i == 1
        legend('Sentman Model', 'IRS Model', 'Location', 'best');
    end
end
sgtitle('Pitch Torque vs Pitch Angle for Individual Faces');

%% Plot Total Torque Results
figure;
plot(rad2deg(pitch_angles__rad), total_torques(:,1,2), 'b-', 'LineWidth', 2);
hold on;
plot(rad2deg(pitch_angles__rad), total_torques(:,2,2), 'r-', 'LineWidth', 2);
grid on;
xlim([0, 360]);
xlabel('Pitch Angle [°]');
ylabel('Total Pitch Torque [Nm]');
legend('Sentman Model', 'IRS Model', 'Location', 'best');
title('Total Pitch Torque vs Pitch Angle');
hold off;

%% Plot Torques by Model (All Bodies)
figure;

subplot(2,1,1);
plot(rad2deg(pitch_angles__rad), squeeze(torques(:,1,2,:)), 'LineWidth', 2);
title('Sentman Model - All Face Torques');
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');
legend(bodies_label, 'Location', 'eastoutside');

subplot(2,1,2);
plot(rad2deg(pitch_angles__rad), squeeze(torques(:,2,2,:)), 'LineWidth', 2);
title('IRS Model - All Face Torques');
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');
legend(bodies_label, 'Location', 'eastoutside');


