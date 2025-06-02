%% LIFT TO DRAG RATIO
% calculate the lift to drag ratio for sentman and the new IRS model
% for a two sided flat plate and plot the results
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

% Import constants from environment:definitions.m
run('environment_definitions.m');

%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'cl_cd_cVAE_A01_flat_and_bird.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% define two sided flat plate
plate = cell(1,1);
plate{1}.vertices_B = zeros(3,3,4);
plate{1}.vertices_B(:,:,1) = [-0.5 0.5 -0.5;
                            -0.5 -0.5 0.5;
                            -1e-10 -1e-10 -1e-10];
plate{1}.vertices_B(:,:,2) = [0.5 0.5 -0.5;
                            0.5 -0.5 0.5;
                            -1e-10 -1e-10 -1e-10];


plate{1}.vertices_B(:,:,3) = [-0.5 0.5 -0.5;
                            -0.5 -0.5 0.5;
                            1e-10 1e-10 1e-10];
plate{1}.vertices_B(:,:,4) = [0.5 0.5 -0.5;
                            0.5 -0.5 0.5;
                            1e-10 1e-10 1e-10]; 



plate{1}.centroids_B = squeeze(mean(plate{1}.vertices_B, 2));
plate{1}.normals_B = [0 0 0 0; 0 0 0 0; -1 -1 1 1];
plate{1}.areas = [0.5 0.5 0.5 0.5];
%rotation hinge point is the centroid of the plate
plate{1}.rotation_hinge_point_B = [0;0;0];
plate{1}.rotation_direction_B = [0;-1;0];
plate{1}.temperatures__K = [300;300;300;300];
plate{1}.energy_accommodation_coefficients = [0.9;0.9;0.9;0.9];

%% calculate aerodynamic forces

%loops
attitude_quarternion_BI = [1;0;0;0];
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        [aerodynamic_force_B__N(:,i,model), ~] = ...
            vleoAerodynamics(...
                attitude_quarternion_BI,...
                rotational_velocity_BI_B__rad_per_s,...
                velocity_I_I__m_per_s,...
                wind_velocity_I_I__m_per_s,...
                density__kg_per_m3,...
                temperature__K,...
                particles_mass__kg,...
                plate,...
                current_angle,...
                temperature_ratio_method,...
                model,...
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
plot(control_surface_angles__rad, lift_to_drag_ratio_sentman, 'DisplayName', 'Sentman Model');
plot(control_surface_angles__rad, lift_to_drag_ratio_new, 'DisplayName', 'New IRS Model');
xlabel('angle of attack [rad]');
ylabel('lift to drag ratio');
title('Lift to Drag Ratio Comparison');
legend;