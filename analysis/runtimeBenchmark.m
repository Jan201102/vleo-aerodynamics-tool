%% RUNTIME BENCHMARK ANALYSIS
% Benchmark runtime performance for different LUT formats:
% Table-based vs polynomial-based lookup tables
import vleo_aerodynamics_core.*
clear;

%% Setup Environment and Geometry
geometry_type = 'shuttlecock';

% Setup for table-based LUT
lut_file_table = 'aerodynamic_coefficients_panel_method.csv';
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data_table] = setup(lut_file_table, geometry_type);

% Setup for polynomial-based LUT  
lut_file_poly = 'aerodynamic_coefficients_panel_method_poly.mat';
[~, ~, ~, ~, ~, lut_data_poly] = setup(lut_file_poly, geometry_type);

%% Runtime Benchmark
num_iterations = 5000;
bodies_rotation_angles__rad = zeros(1, num_bodies);
attitude_quaternion_BI = [1; 0; 0; 0];
model = 2;
execution_times = zeros(2, 1);

fprintf('Starting runtime benchmark with %d iterations...\n', num_iterations);

for m = 1:2
    if m == 1
        fprintf('Testing table-based LUT...\n');
        lut_data = lut_data_table;
    else
        fprintf('Testing polynomial-based LUT...\n');
        lut_data = lut_data_poly;
    end
    
    tic;
    for i = 1:num_iterations
        vleoAerodynamics(...
            attitude_quaternion_BI, ...
            environment_definitions.rotational_velocity_BI_B__rad_per_s, ...
            environment_definitions.velocity_I_I__m_per_s, ...
            environment_definitions.wind_velocity_I_I__m_per_s, ...
            environment_definitions.density__kg_per_m3, ...
            environment_definitions.temperature__K, ...
            environment_definitions.particles_mass__kg, ...
            bodies, ...
            bodies_rotation_angles__rad, ...
            environment_definitions.temperature_ratio_method, ...
            model, ...
            lut_data);
        
        if mod(i, 1000) == 0
            fprintf('Completed %d of %d iterations\n', i, num_iterations);
        end
    end
    
    execution_times(m) = toc;
    fprintf('Execution time for LUT type %d: %.2f seconds\n', m, execution_times(m));
end

%% Display Results
fprintf('\n=== BENCHMARK RESULTS ===\n');
fprintf('Table-based LUT:      %.2f seconds\n', execution_times(1));
fprintf('Polynomial-based LUT: %.2f seconds\n', execution_times(2));

% Calculate relative performance
relative_increase = (execution_times(2) - execution_times(1)) / execution_times(1) * 100;
if relative_increase > 0
    fprintf('Polynomial LUT is %.2f%% SLOWER than table LUT\n', relative_increase);
else
    fprintf('Polynomial LUT is %.2f%% FASTER than table LUT\n', -relative_increase);
end

% Calculate iterations per second
ips_table = num_iterations / execution_times(1);
ips_poly = num_iterations / execution_times(2);
fprintf('Table LUT:      %.0f iterations/second\n', ips_table);
fprintf('Polynomial LUT: %.0f iterations/second\n', ips_poly);

