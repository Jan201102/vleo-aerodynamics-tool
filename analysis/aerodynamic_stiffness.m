%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;
%load environment data
run("environment_definitions.m");
%% Geometry selection parameter
geometry_type = 'shuttlecock'; % Options: 'plate' or 'shuttlecock'
temperature_ratio_method = 1;
%% load lut data
lut_data = load("aerodynamic_coefficients_panel_method_poly.mat");
%% load geometry based on parameter
[bodies, num_bodies, rotation_face_index, x_label] = load_geometry(geometry_type, energy_accommodation, surface_temperature__K);

%% Derivation Configuration
derivation_method = 'central';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_alpha = 1e-4;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% calculate aerodynamic stiffness
num_angles = 100;
control_surface_angles__rad = linspace(0,pi/2, num_angles);
aero_stiffness = zeros(num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        aero_stiffness(i,model) = calculate_aerodynamic_stiffness(derivation_method, ...
                                                                0,...
                                                                delta_alpha, ...
                                                                axis_direction, ...
                                                                torque_component, ...
                                                                bodies, ...
                                                                bodies_rotation_angles__rad, ...
                                                                model, ...
                                                                lut_data);
        fprintf('calcluated point %d of %d\n',i+(model-1)*num_angles,2*num_angles);
    end
end
save("aerodynamic_stiffness_shuttlecock.mat","control_surface_angles__rad","aero_stiffness")
%% plot stiffness of control surface angle for both models in two subplots with a tiled layout
figure;
plot(rad2deg(control_surface_angles__rad), aero_stiffness(:,1), 'b');
hold on;
plot(rad2deg(control_surface_angles__rad), aero_stiffness(:,2), 'r');
grid on;
xlabel(x_label);
ylabel('Aerodynamic Stiffness [Nm/rad]');
legend(sprintf('Sentman \\alpha_E = %.4f',energy_accommodation), 'IRS model',"Location","best");

hold off;

%% Save figure as PNG and EPS
matlab2tikz(sprintf('aerodynamic_stiffness_%s.tex', geometry_type));
%%
figure
plot(rad2deg(control_surface_angles__rad),aero_stiffness(:,1)-aero_stiffness(:,2))
