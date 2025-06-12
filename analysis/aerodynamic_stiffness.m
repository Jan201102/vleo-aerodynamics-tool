%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

%load environment data
run("environment_definitions.m");
%% Geometry selection parameter
geometry_type = 'shuttlecock'; % Options: 'plate' or 'shuttlecock'
energy_accommodation = 1;
%% load geometry based on parameter
if strcmp(geometry_type, 'plate')
    bodies = parametrized_flat_plate(1.0, 1.0, 0.5, 0.0,energy_accommodation);
    showBodies(bodies, [0], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
elseif strcmp(geometry_type, 'shuttlecock')
    bodies = load_from_gmsh(energy_accommodation);
    showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);
    num_bodies = 5;
    rotation_face_index = [2,3,4,5];
else
    error("Invalid geometry_type. Use 'plate' or 'shuttlecock'.");
end

%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'aerodynamic_coefficients_panel_method.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% Derivation Configuration
derivation_method = 'central';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_alpha = 1e-4;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% calculate aerodynamic stiffness
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aero_stiffness = zeros(num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        aero_stiffness(i,model) = calculate_aerodynamic_stiffness(derivation_method, ...
                                                                delta_alpha, ...
                                                                axis_direction, ...
                                                                torque_component, ...
                                                                bodies, ...
                                                                bodies_rotation_angles__rad, ...
                                                                model, ...
                                                                lut);
    end
end

%% plot stiffness of control surface angle for both models in two subplots with a tiled layout
figure;
plot(control_surface_angles__rad, aero_stiffness(:,1), 'b');
hold on;
plot(control_surface_angles__rad, aero_stiffness(:,2), 'r');
grid on;
xlabel('Control Surface Angle [rad]');
ylabel('Aerodynamic Stiffness [Nm/rad]');
legend(sprintf('Sentman \alpha_E = %.2f',energy_accommodation), 'IRS model');
title(sprintf('Aerodynamic Stiffness vs Control Surface Angle for %s geometry', geometry_type));
hold off;

%% Save figure as PNG and EPS
saveas(gcf, sprintf('aerodynamic_stiffness_%s.png', geometry_type));
saveas(gcf, sprintf('aerodynamic_stiffness_%s.eps', geometry_type), 'epsc');
