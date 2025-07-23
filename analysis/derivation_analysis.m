%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

%% load environment data
run("environment_definitions.m");
temperature_ratio_method = 1;
energy_accommodation = 1;
%% load lut data
bodies = load_from_gmsh(energy_accommodation);
lut_data = load_lut("aerodynamic_coefficients_fitted_highres_sentman.csv");

%% Derivation Configuration
derivation_method = 'centrals';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_alpha__rad = 1e-4;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)
alpha_0 = 0;
rotation_angles = [0,pi/4,pi/2];

num_rotation_angles = length(rotation_angles);

%% calculate torque around alpha_0 for both models
num_angles = 101;
attitude_angles = linspace(alpha_0-2*delta_alpha__rad, alpha_0+2*delta_alpha__rad, num_angles);
torque = zeros(num_rotation_angles,num_angles,2,5);
total_torque = zeros(num_rotation_angles,num_angles,2);
for rot_angle = 1:num_rotation_angles
    bodies_rotation_angles__rad = [0;rotation_angles(rot_angle)*ones(4,1)];
    for model = 1:2
        for i = 1:num_angles
            attitude_quaternion_BI = [cos(attitude_angles(i)/2); sin(attitude_angles(i)/2) * axis_direction];
            [~,total_torque_vec,~,torque_vec] = vleoAerodynamics(attitude_quaternion_BI,...
                                                rotational_velocity_BI_B__rad_per_s,...
                                                velocity_I_I__m_per_s,...
                                                wind_velocity_I_I__m_per_s,...
                                                density__kg_per_m3,...
                                                temperature__K,...
                                                particles_mass__kg,...
                                                bodies,...
                                                bodies_rotation_angles__rad,...
                                                temperature_ratio_method,...
                                                model,...
                                                lut_data);
            torque(rot_angle,i,model,:) = torque_vec(torque_component,:);
            total_torque(rot_angle,i,model) = total_torque_vec(torque_component);
            fprintf('calcluated point %d of %d\n',i+(model-1)*num_angles,2*num_angles);
        end
    end
end

%% plot torque curves for both models in a 2x5 grid top row for sentman model, bottom row for new model and each column for a different body
figure;
sgtitle("torque curves for different rotation angles");

% Store handles for legend
legend_handles = [];
legend_labels = {};

for n_body = 1:5
    subplot(2,5,n_body);
    if n_body == 5
        hold on;
        for rot_angle = 1:num_rotation_angles
            h = plot(rad2deg(attitude_angles), squeeze(torque(rot_angle,:,1,5)+torque(rot_angle,:,1,3)), 'DisplayName', sprintf('Rotation Angle %.1f°', rotation_angles(rot_angle)));
            % Store handles only from first iteration
                legend_handles = [legend_handles, h];
                legend_labels{end+1} = sprintf('Rotation Angle %.1f°', rotation_angles(rot_angle));
        end
        title('Body 3+5 (Sentman)');
        xlabel('Attitude Angle (°)');
        ylabel('Torque (Nm)');
        
        subplot(2,5,n_body+5);
        hold on;
        for rot_angle = 1:num_rotation_angles
            plot(rad2deg(attitude_angles), squeeze(torque(rot_angle,:,2,5)+torque(rot_angle,:,2,3)), 'DisplayName', sprintf('Rotation Angle %.1f°', rotation_angles(rot_angle)));
        end
        title('Body 3+5 (New Model)');
        xlabel('Attitude Angle (°)');
        ylabel('Torque (Nm)');
    else
        hold on;
        for rot_angle = 1:num_rotation_angles
            plot(rad2deg(attitude_angles), squeeze(torque(rot_angle,:,1,n_body)), 'DisplayName', sprintf('Rotation Angle %.1f°', rotation_angles(rot_angle)));
        end
        title(sprintf('Body %d (Sentman)', n_body));
        xlabel('Attitude Angle (°)');
        ylabel('Torque (Nm)');
        
        subplot(2,5,n_body+5);
        hold on;
        for rot_angle = 1:num_rotation_angles
            plot(rad2deg(attitude_angles), squeeze(torque(rot_angle,:,2,n_body)), 'DisplayName', sprintf('Rotation Angle %.1f°', rotation_angles(rot_angle)));
        end
        title(sprintf('Body %d (New Model)', n_body));
        xlabel('Attitude Angle (°)');
        ylabel('Torque (Nm)');
    end
end

% Create shared legend positioned on the right side
lgd = legend(legend_handles, legend_labels, 'Position', [0.92, 0.3, 0.06, 0.4]);
lgd.Title.String = 'Rotation Angles';

%% calculate torque at derivation points
torque_at_nodes = zeros(num_rotation_angles,2,2);
nodes = [alpha_0-delta_alpha__rad, alpha_0+delta_alpha__rad];
for rot_angle = 1:num_rotation_angles
    bodies_rotation_angles__rad = [0;rotation_angles(rot_angle)*ones(4,1)];
    for model = 1:2
        for i = 1:2
            alpha = nodes(i);
            attitude_quaternion_BI = [cos(alpha/2); sin(alpha/2) * axis_direction];
            [~,torque_vec,~,~] = vleoAerodynamics(attitude_quaternion_BI,...
                                                rotational_velocity_BI_B__rad_per_s,...
                                                velocity_I_I__m_per_s,...
                                                wind_velocity_I_I__m_per_s,...
                                                density__kg_per_m3,...
                                                temperature__K,...
                                                particles_mass__kg,...
                                                bodies,...
                                                bodies_rotation_angles__rad,...
                                                temperature_ratio_method,...
                                                model,...
                                                lut_data);
            torque_at_nodes(rot_angle,i,model) = torque_vec(torque_component);
        end
    end
end
derivatives = (torque_at_nodes(:,2,:)-torque_at_nodes(:,1,:))/(2*delta_alpha__rad);
disp(derivatives);
%% plot torques and torques at derivation nodes, each model has a subplot
figure;
sgtitle("torques and derivation points");
subplot(1,2,1);
hold on;
for rot_angle = 1:num_rotation_angles
    if rot_angle == 1
        plot(rad2deg(attitude_angles), total_torque(rot_angle,:,1), 'DisplayName', sprintf('Rotation Angle %.1f', rotation_angles(rot_angle)));
        plot(rad2deg(nodes), squeeze(torque_at_nodes(rot_angle,:,1)), 'ro', 'DisplayName', 'Derivation Points');
    else
        plot(rad2deg(attitude_angles), total_torque(rot_angle,:,1), 'DisplayName', sprintf('Rotation Angle %.1f', rotation_angles(rot_angle)));
        plot(rad2deg(nodes), squeeze(torque_at_nodes(rot_angle,:,1)), 'ro', 'HandleVisibility', 'off');
    end
end
title('sentman Model');
xlabel('Attitude Angle (°)');
ylabel('Torque (Nm)');
legend;
subplot(1,2,2);
hold on;
for rot_angle = 1:num_rotation_angles
    if rot_angle == 1
        plot(rad2deg(attitude_angles), total_torque(rot_angle,:,2), 'DisplayName', sprintf('Rotation Angle %.1f', rotation_angles(rot_angle)));
        plot(rad2deg(nodes), squeeze(torque_at_nodes(rot_angle,:,2)), 'ro', 'DisplayName', 'Derivation Points');
    else
        plot(rad2deg(attitude_angles), total_torque(rot_angle,:,2), 'DisplayName', sprintf('Rotation Angle %.1f', rotation_angles(rot_angle)));
        plot(rad2deg(nodes), squeeze(torque_at_nodes(rot_angle,:,2)), 'ro', 'HandleVisibility', 'off');
    end
end
title('new Model');   
xlabel('Attitude Angle (°)');
ylabel('Torque (Nm)');
legend;