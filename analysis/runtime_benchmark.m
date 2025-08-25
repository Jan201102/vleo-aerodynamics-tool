%% RUNTIME BENCHMARK
% benchmark runtime for different LUT formats
import vleo_aerodynamics_core.*
clear;

%% Setup environment and geometry
geometry_type = 'shuttlecock';
% Setup for table-based LUT
lut_file_table = "aerodynamic_coefficients_panel_method.csv";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data_table] = setup(lut_file_table, geometry_type);

% Setup for polynomial-based LUT  
lut_file_poly = "aerodynamic_coefficients_panel_method_poly.mat";
[~, ~, ~, ~, ~, lut_data_poly] = setup(lut_file_poly, geometry_type);

%% calculate aerodynamic stiffness
num_iterations = 5000;
bodies_rotation_angles__rad = zeros(1, num_bodies);
attitude_quaternion_BI = [1;0;0;0];
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
            attitude_quaternion_BI,...
            environment_definitions.rotational_velocity_BI_B__rad_per_s,...
            environment_definitions.velocity_I_I__m_per_s,...
            environment_definitions.wind_velocity_I_I__m_per_s,...
            environment_definitions.density__kg_per_m3,...
            environment_definitions.temperature__K,...
            environment_definitions.particles_mass__kg,...
            bodies,...
            bodies_rotation_angles__rad,...
            environment_definitions.temperature_ratio_method,...
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

