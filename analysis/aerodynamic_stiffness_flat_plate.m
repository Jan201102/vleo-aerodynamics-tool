%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

%load environment data
run("environment_definitions.m");

%% flat plate
% Create a flat plate using the parametrized function
% Parameters: x_dim=1, y_dim=1, cog_x=0, cog_y=0 (unit square centered at origin)
bodies = parametrized_flat_plate(1.0, 1.0, -0.5, 0.0);

%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'cl_cd_cVAE_A01_flat_and_bird.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% Derivation Configuration
derivation_method = 'seven_point_stencil';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_alpha = pi/1800000;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% show bodies
showBodies(bodies, [0], 0.75, 0.25);

%% calculate aerodynamic stiffness
aero_stiffness = zeros(1,3);
for i = 1:3
    aero_stiffness(i) = calculate_aerodynamic_stiffness('seven_point_stencil', ...
                                                    pi/1000, ...
                                                    [0; 1; 0], ...
                                                    2, ...
                                                    bodies, ...
                                                    [0], ...
                                                    i, ...
                                                    lut);
end
disp('Aerodynamic Stiffness:');
disp(aero_stiffness);