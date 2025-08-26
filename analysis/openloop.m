%% OPEN LOOP ANALYSIS
% Calculate the aerodynamic stiffness and damping for open loop analysis
% Analyze natural frequencies and damping ratios for different control surface angles
import vleo_aerodynamics_core.*
clear;

%% Setup Environment, Geometry, and LUT
geometry_type = 'shuttlecock';
lut_file = 'aerodynamic_coefficients_panel_method_poly.mat';
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);
I_yy = 0.0375; % Moment of inertia around y-axis

%% Derivation Configuration
derivation_method = 'central'; 
delta_q_rad_per_s = 1e2;          % Step size for numerical differentiation [rad/s]
delta_alpha = 1e-4;               % Step size for numerical differentiation [rad]
axis_direction = [0; 1; 0];       % Rotation axis (y-axis for pitch)
torque_component = 2;             % Component index (2 = y-component for pitch)

%% Calculate Aerodynamic Stiffness and Damping
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
ol_parameters = zeros(num_angles, 2, 2);  % [angles, models, parameters]

for model = 1:2
    fprintf('Processing model %d...\n', model);
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        % Calculate aerodynamic damping
        C_M_da = calculate_aerodynamic_damping(...
            derivation_method, ...
            delta_q_rad_per_s, ...
                                                axis_direction, ...
            torque_component, ...
            bodies, ...
            bodies_rotation_angles__rad, ...
            model, ...
            lut_data, ...
            environment_definitions);
        
        % Calculate aerodynamic stiffness
        C_M_a = calculate_aerodynamic_stiffness(...
            derivation_method, ...
            0, ...
            delta_alpha, ...
            axis_direction, ...
            torque_component, ...
            bodies, ...
            bodies_rotation_angles__rad, ...
            model, ...
            lut_data, ...
            environment_definitions);

        % Calculate open-loop parameters
        omega_0 = sqrt(-C_M_a / I_yy);                    % Natural frequency [rad/s]
        zeta = -C_M_da / (2 * I_yy * omega_0);           % Damping ratio [-]
        
        ol_parameters(i, model, 1) = omega_0;
        ol_parameters(i, model, 2) = zeta;
        
        if mod(i, 10) == 0
            fprintf('Calculated point %d of %d\n', i + (model-1)*num_angles, 2*num_angles);
        end
    end
end

%% Save Calculated Open-Loop Parameters
save('ol_parameters.mat', 'ol_parameters', 'control_surface_angles__rad');

%% Plot Natural Frequency and Damping Ratio
figure;
t = tiledlayout(2, 1);

% Natural frequency plot
nexttile;
plot(rad2deg(control_surface_angles__rad), squeeze(ol_parameters(:,1,1)), 'b-');
hold on;
plot(rad2deg(control_surface_angles__rad), squeeze(ol_parameters(:,2,1)), 'r-');
grid on;
xlabel(x_label);
ylabel('Natural Frequency \omega_0 [rad/s]');
title('Natural Frequency for Different Models');
legend('Sentman', 'IRS', 'Location', 'southeast');

% Damping ratio plot
nexttile;
plot(rad2deg(control_surface_angles__rad), squeeze(ol_parameters(:,1,2)), 'b-');
hold on;
plot(rad2deg(control_surface_angles__rad), squeeze(ol_parameters(:,2,2)), 'r-');
grid on;
xlabel(x_label);
ylabel('Damping Ratio \zeta [-]');
title('Damping Ratio for Different Models');
legend('Sentman', 'IRS', 'Location', 'southeast');

% Save plot
matlab2tikz('ol_parameters_shuttlecock.tex');

%% Compare Natural and Damped Frequencies
% Calculate damped frequency and relative difference
omega_d = ol_parameters(:,:,1) .* sqrt(1 - ol_parameters(:,:,2).^2);
omega_0 = ol_parameters(:,:,1);
relative_difference = (omega_0 - omega_d) ./ omega_0;

figure;
plot(rad2deg(control_surface_angles__rad), relative_difference(:,1), 'b-');
hold on;
plot(rad2deg(control_surface_angles__rad), relative_difference(:,2), 'r-');
grid on;
xlabel(x_label);
ylabel('Relative Difference [-]');
title('Relative Difference between \omega_0 and \omega_d');
legend('Sentman', 'IRS', 'Location', 'southeast');


 
