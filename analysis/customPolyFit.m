%% CUSTOM POLYNOMIAL FIT ANALYSIS
% Process aerodynamic coefficient data and create polynomial fits
% for lift and drag coefficients with proper constraints
clear; clc; close all;

%% Load and Process Data
csv_file_in = 'aerodynamic_coefficients_panel_method.csv';
data = readtable(csv_file_in);

% Extract data
AoA = data.AoA;
C_l_ram = data.C_l_ram;

% Ensure 90° lift value is 0 (physical constraint)
C_l_ram(end) = 0;
C_d_ram = data.C_d_ram;
C_l_wake = data.C_l_wake;
C_d_wake = data.C_d_wake;

%% Combine Ram and Wake Data
% Reorder and concatenate data to cover full angle range
AoA_total = [flip(-AoA(2:end)); AoA];
C_l_total = [flip(-C_l_wake(2:end)); C_l_ram];
C_d_total = [flip(C_d_wake(2:end)); C_d_ram];
data = [AoA_total, C_l_total, C_d_total];

%% Plot Raw Data
figure;
hold on;
scatter(AoA_total, C_l_total, '.', 'DisplayName', 'C_l Data', 'SizeData', 50);
scatter(AoA_total, C_d_total, '.', 'DisplayName', 'C_d Data', 'SizeData', 50);
xlabel('Angle of Attack [°]');
ylabel('Aerodynamic Coefficients [-]');
title('Raw Aerodynamic Coefficient Data');
legend('Location', 'best');
grid on;

%% Polynomial Degree Analysis
% Define boundaries for piecewise polynomial segments
boundaries = [-24, 0];  % [boundary_1_2, boundary_2_3]

% Analyze polynomial degrees from 4 to 20 for segments 2 and 3
degree_range = 4:20;
mse_cl_values = zeros(size(degree_range));
mse_cd_values = zeros(size(degree_range));

% Get segment 2 and 3 indices for evaluation
segment_2_3_idx = (AoA_total > boundaries(1) & AoA_total <= boundaries(2)) | ...
                  (AoA_total > boundaries(2));
n_points_eval = sum(segment_2_3_idx);  % Number of points for BIC calculation

fprintf('Evaluating polynomial degrees from %d to %d...\n', min(degree_range), max(degree_range));

for i = 1:length(degree_range)
    degree = degree_range(i);
    
    % Polynomial degrees for each segment: [constant, degree, degree]
    poly_degrees = [0, degree, degree];
    
    % Constraint points for C_l (lift coefficient must be 0 at 90°)
    constraint_points_Cl = [90, 0];
    
    % Fit polynomials
    global_piecewise_coeff_cl = piecewise_poly_fit_with_constraints(AoA_total, C_l_total, boundaries, poly_degrees, constraint_points_Cl, 1);
    global_piecewise_coeff_cd = piecewise_poly_fit_with_constraints(AoA_total, C_d_total, boundaries, poly_degrees, [], 1);
    
    % Convert to piecewise polynomials
    pp_Cl = coeffs_to_mkpp(global_piecewise_coeff_cl, boundaries, poly_degrees);
    pp_Cd = coeffs_to_mkpp(global_piecewise_coeff_cd, boundaries, poly_degrees);
    
    % Calculate predictions for segments 2 and 3 only
    Cl_pred = ppval(pp_Cl, AoA_total(segment_2_3_idx));
    Cd_pred = ppval(pp_Cd, AoA_total(segment_2_3_idx));
    
    % Calculate MSE for segments 2 and 3 only
    mse_cl_values(i) = mean((C_l_total(segment_2_3_idx) - Cl_pred).^2);
    mse_cd_values(i) = mean((C_d_total(segment_2_3_idx) - Cd_pred).^2);

    fprintf('Degree %d: MSE_Cl = %.6f, MSE_Cd = %.6f\n', ...
        degree, mse_cl_values(i), mse_cd_values(i));
end

%% Plot MSE vs polynomial degree
figure;
hold on;
plot(degree_range, mse_cl_values);
plot(degree_range, mse_cd_values);
xlabel('Polynomial Degree');
ylabel('MSE');
title('MSE vs Polynomial Degree');
legend('C_l','C_d')
grid on;
set(gca, 'YScale', 'log'); % Use log scale for better visualization

% Find optimal degrees based on MSE
[min_mse_cl, min_idx_cl] = min(mse_cl_values);
[min_mse_cd, min_idx_cd] = min(mse_cd_values);
optimal_degree_cl = degree_range(min_idx_cl);
optimal_degree_cd = degree_range(min_idx_cd);

fprintf('\nOptimal degrees (MSE-based):\n');
fprintf('C_l: Degree %d (MSE = %.6f)\n', optimal_degree_cl, min_mse_cl);
fprintf('C_d: Degree %d (MSE = %.6f)\n', optimal_degree_cd, min_mse_cd);

%% Fit with optimal degrees and visualize
poly_degrees_c_l = [0, 8,8];
poly_degrees_c_d = [0,9,9];
constraint_points_Cl = [90,0];

global_piecewise_coeff_cl = piecewise_poly_fit_with_constraints(AoA_total, C_l_total, boundaries, poly_degrees_c_l, constraint_points_Cl, 1);
global_piecewise_coeff_cd = piecewise_poly_fit_with_constraints(AoA_total, C_d_total, boundaries, poly_degrees_c_d, [], 1);
pp_Cl = coeffs_to_mkpp(global_piecewise_coeff_cl, boundaries, poly_degrees_c_l);
pp_Cd = coeffs_to_mkpp(global_piecewise_coeff_cd, boundaries, poly_degrees_c_d);
save('aerodynamic_coefficients_panel_method_poly.mat','pp_Cd','pp_Cl');
%% Test and plot the results
AoA_test = linspace(-90, 90, 1000);
Cl_combined = ppval(pp_Cl, AoA_test);
Cd_combined = ppval(pp_Cd, AoA_test);

% Plot C_l
figure;
subplot(2,1,1);
hold on;
scatter(AoA_total, C_l_total, 'DisplayName', 'Original C_l Data');
plot(AoA_test, Cl_combined, 'DisplayName', 'Combined Piecewise Polynomial');
xlabel('Angle of Attack (degrees)');
ylabel('C_l');
title('C_l Piecewise Polynomial Fit');
legend('Location', 'best');
grid on;

% Plot C_d
subplot(2,1,2);
hold on;
scatter(AoA_total, C_d_total, 'DisplayName', 'Original C_d Data');
plot(AoA_test, Cd_combined, 'DisplayName', 'Combined Piecewise Polynomial');

xlabel('Angle of Attack (degrees)');
ylabel('C_d');
title('C_d Piecewise Polynomial Fit');
legend('Location', 'best');
grid on;

%% Plot residuals for both coefficients
figure;
subplot(1,2,1);
C_l_total_pred = ppval(pp_Cl, AoA_total);
scatter(AoA_total, C_l_total_pred - C_l_total);
grid on;
xlabel('Angle of Attack (degrees)');
ylabel('Residuals (C_l)');

subplot(1,2,2);
C_d_total_pred = ppval(pp_Cd, AoA_total);
scatter(AoA_total, C_d_total_pred - C_d_total);
grid on;
xlabel('Angle of Attack (degrees)');
ylabel('Residuals (C_d)');