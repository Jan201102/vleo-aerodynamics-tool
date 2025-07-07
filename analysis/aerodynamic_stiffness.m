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
lut_data = load_lut("aerodynamic_coefficients_panel_method.csv");


%% load geometry based on parameter
if strcmp(geometry_type, 'plate')
    bodies = parametrized_flat_plate(1, 1, 0.5, 0.0,true,energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [°]";
elseif strcmp(geometry_type, 'shuttlecock')
    bodies = load_from_gmsh(energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);
    num_bodies = 5;
    rotation_face_index = [2,3,4,5];
    x_label = "control surface angle [°]";
elseif strcmp(geometry_type, 'shuttlecock_wing')
    bodies = load_shuttlecock_wing(energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [°]";
elseif strcmp(geometry_type, 'shuttlecock_wing_new')
    bodies_all = load_from_gmsh(energy_accommodation,surface_temperature__K);
    bodies = cell(1,1);
    bodies{1} = bodies_all{3};
    showBodies(bodies, [pi/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [°]";
else
    error("Invalid geometry_type. Use 'plate' or 'shuttlecock'.");
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
                                                                lut_data);
        fprintf('calcluated point %d of %d\n',i+(model-1)*num_angles,2*num_angles);
    end
end

%% plot stiffness of control surface angle for both models in two subplots with a tiled layout
figure;
plot(rad2deg(control_surface_angles__rad), aero_stiffness(:,1), 'b');
hold on;
plot(rad2deg(control_surface_angles__rad), aero_stiffness(:,2), 'r');
grid on;
xlabel(x_label);
ylabel('Aerodynamic Stiffness [Nm/rad]');
legend(sprintf('Sentman \\alpha_E = %.4f',energy_accommodation), 'IRS model');

hold off;

%% Save figure as PNG and EPS
matlab2tikz(sprintf('aerodynamic_stiffness_%s.tex', geometry_type));
