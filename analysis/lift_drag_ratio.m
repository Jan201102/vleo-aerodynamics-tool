%% LIFT TO DRAG RATIO
% calculate the lift to drag ratio for sentman and the new IRS model
% for a two sided flat plate and plot the results
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

% Import constants from environment:definitions.m
run('environment_definitions.m');

%% Geometry selection parameter
% Set to 'plate' for flat plate geometry or 'shuttlecock' for shuttlecock model
geometry_type = 'shuttlecock'; % Options: 'plate' or 'shuttlecock'
energy_accomodation = 0.9;
%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'aerodynamic_coefficients_panel_method.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% load geometry based on parameter
if strcmp(geometry_type, 'plate')
    bodies = parametrized_flat_plate(0.13333, 0.098, 0.0, 0.0,true,energy_accomodation);
    showBodies(bodies, [0], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
elseif strcmp(geometry_type, 'shuttlecock')
    bodies = load_from_gmsh(energy_accomodation);
    showBodies(bodies, [0,0/4,0/4,pi/4,0/4], 0.75, 0.25);
    num_bodies = 5;
    rotation_face_index = 3;
elseif strcmp(geometry_type, 'shuttlecock_wing')
    bodies = load_shuttlecock_wing();
    showBodies(bodies, [0/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
elseif strcmp(geometry_type, 'shuttlecock_wing_new')
    bodies_all = load_from_gmsh(energy_accomodation);
    bodies = cell(1,1);
    bodies{1} = bodies_all{3};
    showBodies(bodies, [pi/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
else
    error("Invalid geometry_type. Use 'plate' or 'shuttlecock'.");
end

%% calculate aerodynamic forces

%loops
attitude_quarternion_BI = [1;0;0;0];
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        [aerodynamic_force_B__N(:,i,model), ~] = ...
            vleoAerodynamics(...
                attitude_quarternion_BI,...
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
                1,...
                lut);
    end
end

aerodynamic_force_sentman__N = aerodynamic_force_B__N(:,:,1);
aerodynamic_force_new__N = aerodynamic_force_B__N(:,:,2);
% calculate lift to drag ratio
lift_sentman__N = aerodynamic_force_sentman__N(3,:);
drag_sentman__N = -aerodynamic_force_sentman__N(1,:);
lift_new__N = aerodynamic_force_new__N(3,:);
drag_new__N = -aerodynamic_force_new__N(1,:);
lift_to_drag_ratio_sentman = lift_sentman__N./drag_sentman__N;
lift_to_drag_ratio_new = lift_new__N./drag_new__N;

%% plot results
figure;
hold on;
grid on;
plot(control_surface_angles__rad, lift_to_drag_ratio_sentman,"b", 'DisplayName', sprintf('Sentman Model \\alpha_E = %.4f',energy_accomodation));
plot(control_surface_angles__rad, lift_to_drag_ratio_new,"r", 'DisplayName', 'New IRS Model');
xlabel('control surface angle [rad]');
ylabel('lift to drag ratio');
title(sprintf('Lift to Drag Ratio Comparison for %s geometry',geometry_type));
legend;

%% Save figure as PNG and EPS
saveas(gcf, sprintf('lift_drag_ratio_%s.png', geometry_type));
saveas(gcf, sprintf('lift_drag_ratio_%s.eps', geometry_type), 'epsc');

%% plot force envelopes for both models
figure;
hold on;
grid on;
plot3(aerodynamic_force_B__N(1,:,1), aerodynamic_force_B__N(2,:,1), aerodynamic_force_B__N(3,:,1), 'b', 'DisplayName', sprintf('Sentman Model \\alpha_E = %.4f',energy_accomodation));
plot3(aerodynamic_force_B__N(1,:,2), aerodynamic_force_B__N(2,:,2), aerodynamic_force_B__N(3,:,2), 'r', 'DisplayName', 'New IRS Model');
xlabel('X Force [N]');
ylabel('Y Force [N]');
zlabel('Z Force [N]');
title(sprintf('Aerodynamic Force Envelopes for %s geometry', geometry_type));
view([0 -1 0])
legend;