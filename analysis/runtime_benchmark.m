%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;
run('environment_definitions.m');
%% load lut data
lut_data_table = load_lut("aerodynamic_coefficients_panel_method.csv");
lut_data_poly = load("aerodynamic_coefficients_panel_method_poly.mat");

bodies = load_from_gmsh();
num_bodies = 5;

%% calculate aerodynamic stiffness
num_iterations = 100;
bodies_rotation_angles__rad = zeros(1, num_bodies);
attitude_quarternion_BI = [1;0;0;0];
model = 2;
execution_times = zeros(2,1);
for m = 1:2
    tic;
    for i = 1:num_iterations
        if m == 1
            lut_data = lut_data_table;
        else
            lut_data = lut_data_poly;
        end
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
            lut_data);
        fprintf('calcluated point %d of %d\n',i+(m-1)*num_iterations,2*num_iterations);
    end
    execution_times(m) = toc;
    fprintf('Execution time for model %d: %.2f seconds\n', m, execution_times(m));
end
%relative increase in execution time
relative_increase = (execution_times(2) - execution_times(1)) / execution_times(1) * 100;
fprintf('Relative increase in execution time: %.2f%%\n', relative_increase);

