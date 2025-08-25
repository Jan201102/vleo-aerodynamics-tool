% Aerodynamics
import vleo_aerodynamics_core.*;

%% Setup environment, geometry, and LUT
geometry_type = 'unit_cube';
lut_file = "aerodynamic_coefficients_panel_method_poly.mat";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);
bodies_label = ["bottom","left","top","right","back","front"];

%% Calculate torque curves

% Set control surface angles to 0 for shuttlecock geometry
bodies_rotation_angles = zeros(num_bodies);
num_pitch_angles = 201;
pitch_angles__rad = linspace(0,2*pi, num_pitch_angles);
%pitch_angles__rad = [pi/8];
torques = zeros(num_pitch_angles, 2,3,num_bodies); % 2 models
total_torques = zeros(num_pitch_angles, 2,3); % 2 models
axis_direction = [0; 1; 0]; % Rotation axis (y-axis for pitch)
% Calculate torque for each pitch angle and model (only shuttlecock)
for model = 1:2
    for i = 1:num_pitch_angles
        pitch_angle = pitch_angles__rad(i);
        %bodies_rotation_angles = [pitch_angle];
        % Calculate torque for shuttlecock geometry
        attitude_quaternion_BI = [cos(pitch_angle/2); sin(pitch_angle/2) * axis_direction];
        %attitude_quaternion_BI = [1;0;0;0];
        [~, total_torque_vec,~,torque_vec] = vleoAerodynamics(...
            attitude_quaternion_BI,...
            environment_definitions.rotational_velocity_BI_B__rad_per_s,...
            environment_definitions.velocity_I_I__m_per_s,...
            environment_definitions.wind_velocity_I_I__m_per_s,...
            environment_definitions.density__kg_per_m3,...
            environment_definitions.temperature__K,...
            environment_definitions.particles_mass__kg,...
            bodies,...
            bodies_rotation_angles,...
            environment_definitions.temperature_ratio_method,...
            model,...
            lut_data);
        torques(i, model,:,:) = torque_vec;
        total_torques(i, model,:) = total_torque_vec';
    end
end
%% plot pitch torques for each body in a 3x2 grid
figure;
for i = 1:num_bodies
    subplot(3,2,i);
    hold on;
    % Plot torques for each model
    plot(rad2deg(pitch_angles__rad), torques(:,1,2,i), 'b-', 'LineWidth', 2);
    plot(rad2deg(pitch_angles__rad), torques(:,2,2,i), 'r-', 'LineWidth', 2);
    
    grid on;
    xlabel('global Pitch Angle [°]');
    ylabel('local Pitch Torque [Nm]');
    title(bodies_label(i));
end
%% Plot results
c_total_torques = squeeze(sum(torques, 4)); % Sum over all bodies
figure;
plot(rad2deg(pitch_angles__rad), total_torques(:,1,2), 'b-');
hold on;
plot(rad2deg(pitch_angles__rad), total_torques(:,2,2), 'r-');
grid on
xlim([0 360])
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');
legend("Sentman","IRS")
title('Pitch Torque vs Pitch Angle ');
hold off;
%matlab2tikz('torque_curves_unit_cube.tex')
%% plot pitch torques for each model of all bodies in a 2x1 grid
figure;
subplot(2,1,1);
plot(rad2deg(pitch_angles__rad), squeeze(torques(:,1,2,:)), 'b-', 'LineWidth', 2);
title("sentman")
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');
subplot(2,1,2);
plot(rad2deg(pitch_angles__rad), squeeze(torques(:,2,2,:)), 'r-', 'LineWidth', 2);
title("new model")
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');


