import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

% Import constants from environment:definitions.m
run('environment_definitions.m');

%% Geometry selection parameter
geometry_type = 'shuttlecock'; % Options: 'plate' or 'shuttlecock'
temperature_ratio_method = 1;

%% load lut data
lut_data = load_lut("aerodynamic_coefficients_fitted_highres.csv");

%% load geometry based on parameter
if strcmp(geometry_type, 'plate')
    bodies = parametrized_flat_plate(1, 1, 0.5, 0.0,true,energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [deg]";
elseif strcmp(geometry_type, 'shuttlecock')
    bodies = load_from_gmsh(energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);
    num_bodies = 5;
    rotation_face_index = [2,3,4,5];
    x_label = "control surface angle [deg]";
elseif strcmp(geometry_type, 'shuttlecock_wing')
    bodies = load_shuttlecock_wing(energy_accommodation,surface_temperature__K);
    showBodies(bodies, [0/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [deg]";
elseif strcmp(geometry_type, 'shuttlecock_wing_new')
    bodies_all = load_from_gmsh(energy_accommodation,surface_temperature__K);
    bodies = cell(1,1);
    bodies{1} = bodies_all{3};
    showBodies(bodies, [pi/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
    x_label = "angle of attack [deg]";
else
    error("Invalid geometry_type. Use 'plate' or 'shuttlecock'.");
end

%% Define delta_alpha values to analyze
delta_alphas = [1e-4, 1e-3, 1e-2,1e-1]; % radians (reduced range for angle sweep)
num_deltas = length(delta_alphas);

% Define test angle - only evaluate at 0 degrees for initial analysis
test_angle__rad = 0; % 0 degrees

% Define control surface angles for stiffness vs angle analysis
control_surface_angles__rad = linspace(0, pi/2, 21); % 0 to 90 degrees
num_control_angles = length(control_surface_angles__rad);

% Define numerical differentiation methods to compare
methods = {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'};
num_methods = length(methods);

% Storage for results
stiffness_vs_angle = nan(num_control_angles, num_deltas, num_methods); % For angle sweep analysis

axis_dir_B = [0; 1; 0]; % Rotation around y-axis (pitch)
component = 2; % y-component of torque (pitch torque)
model = 2; % Use IRS model

%% Calculate stiffness vs control surface angle for different delta_alpha values
fprintf('\nAnalyzing stiffness vs control surface angle for different delta_alpha values...\n');

for delta_idx = 1:num_deltas
    delta = delta_alphas(delta_idx);
    fprintf('Computing stiffness vs angle for delta_alpha = %.2e rad\n', delta);
    
    for angle_idx = 1:num_control_angles
        current_angle = control_surface_angles__rad(angle_idx);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        for method_idx = 1:num_methods
            method = methods{method_idx};
            
            try
                stiffness = calculate_aerodynamic_stiffness(...
                    method, delta, axis_dir_B, component, ...
                    bodies, bodies_rotation_angles__rad, model, lut_data);
                stiffness_vs_angle(angle_idx, delta_idx, method_idx) = stiffness;
            catch ME
                % Silently handle errors for cleaner output
                disp(ME.message)
                stiffness_vs_angle(angle_idx, delta_idx, method_idx) = NaN;
            end
        end
    end
end

%% Calculate absolute changes in stiffness with respect to delta_alpha
% Calculate absolute differences between consecutive delta_alpha values
abs_changes = nan(num_control_angles, num_deltas-1, num_methods);

for method_idx = 1:num_methods
    for angle_idx = 1:num_control_angles
        for delta_idx = 1:num_deltas-1
            stiffness_current = stiffness_vs_angle(angle_idx, delta_idx, method_idx);
            stiffness_next = stiffness_vs_angle(angle_idx, delta_idx+1, method_idx);
            
            if ~isnan(stiffness_current) && ~isnan(stiffness_next)
                abs_changes(angle_idx, delta_idx, method_idx) = ...
                    abs(stiffness_next - stiffness_current);
            end
        end
    end
end

% Calculate average absolute changes across control angles
avg_abs_changes = squeeze(mean(abs_changes, 1, 'omitnan'));

%% Plot average absolute changes across all methods
figure();
hold on;
grid on;

method_colors = lines(num_methods);
delta_pairs = delta_alphas(1:end-1);

for method_idx = 1:num_methods
    plot(delta_pairs, avg_abs_changes(:, method_idx), 'o-', ...
             'Color', method_colors(method_idx,:), 'LineWidth', 2, ...
             'DisplayName', methods{method_idx});
end

xlabel('Delta Alpha [rad]');
ylabel('Average Absolute Change [Nm/rad]');
title('Average Absolute Changes in Stiffness vs Delta Alpha');
legend('Location', 'best');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
%% plot stiffness vs control surface angle for each method in a tiled layout
figure();
tl = tiledlayout(num_methods, 1, 'TileSpacing', 'compact');
for method_idx = 1:num_methods
    nexttile;
    hold on;
    grid on;
    
    for delta_idx = 1:num_deltas
        plot(control_surface_angles__rad, ...
             squeeze(stiffness_vs_angle(:, delta_idx, method_idx)), ...
             'DisplayName', sprintf('Delta = %.2e rad', delta_alphas(delta_idx)));
    end
    
    title(sprintf('Stiffness vs Control Surface Angle (%s)', methods{method_idx}));
    xlabel(x_label);
    ylabel('Aerodynamic Stiffness [Nm/rad]');
    legend('Location', 'best');
end