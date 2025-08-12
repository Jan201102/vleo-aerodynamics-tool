% Aerodynamics
import vleo_aerodynamics_core.*;
%load environment data
run("environment_definitions.m");

%load bodies
%bodies = load_box(...
%    energy_accommodation,...
%    surface_temperature__K...
%);

% bodies = parametrized_flat_plate(...
%     1, 1,...
%     [0.5;0;0],...
%     false,energy_accommodation,surface_temperature__K);
 bodies = load_box_truncated(...
     energy_accommodation,...
     surface_temperature__K...
 );
 bodies_label = ["bottom","left","top","right","back","front"]
num_bodies = length(bodies);

%load lut data
lut_data = load("aerodynamic_coefficients_panel_method_poly.mat");
%%

showBodies(bodies,[0,0,0,0,0,0])

%% Calculate torque curves

% Set control surface angles to 0 for shuttlecock geometry
bodies_rotation_angles = pi/2 * zeros(num_bodies);
num_pitch_angles = 101;
pitch_angles__rad = linspace(-pi/4, pi/4, num_pitch_angles);
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
            rotational_velocity_BI_B__rad_per_s,...
            velocity_I_I__m_per_s,...
            wind_velocity_I_I__m_per_s,...
            density__kg_per_m3,...
            temperature__K,...
            particles_mass__kg,...
            bodies,...
            bodies_rotation_angles,...
            temperature_ratio_method,...
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
subplot(2,1,1);
plot(rad2deg(pitch_angles__rad), total_torques(:,1,2), 'b-', 'LineWidth', 2);
hold on;
plot(rad2deg(pitch_angles__rad), c_total_torques(:,1,2), 'b--', 'LineWidth', 2);
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');

subplot(2,1,2)
plot(rad2deg(pitch_angles__rad), total_torques(:,2,2), 'r-', 'LineWidth', 2);
hold on;
plot(rad2deg(pitch_angles__rad), c_total_torques(:,2,2), 'r--', 'LineWidth', 2);
grid on;
xlabel('Pitch Angle [°]');
ylabel('Pitch Torque [Nm]');

title('Pitch Torque vs Pitch Angle ');
hold off;

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


