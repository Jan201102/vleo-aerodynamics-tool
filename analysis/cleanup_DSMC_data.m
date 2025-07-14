%% Aerodynamic Coefficient Function Fitting with Smart Wake Handling
clear; clc; close all;
%% Load CSV data
csv_file_in = 'aerodynamic_coefficients_panel_method.csv';
output_filename = 'aerodynamic_coefficients_panel_method_clean.csv';
data = readtable(csv_file_in);

% Extract data
AoA = data.AoA;
C_l_ram = data.C_l_ram;
C_d_ram = data.C_d_ram;
C_l_wake = data.C_l_wake;
C_d_wake = data.C_d_wake;
%% and reorder concatenate data
AoA_total = [flip(-AoA(2:end));AoA];
C_l_total = [flip(C_l_wake(2:end));C_l_ram];
C_d_total = [flip(C_d_wake(2:end));C_d_ram];

%% Analyze wake data to find cutoff point
% Find where wake coefficients become effectively zero
wake_threshold = 5e-3;
C_l_nonzero_idx = find(abs(C_l_total) > wake_threshold, 1, 'first');
C_d_nonzero_idx = find(abs(C_d_total) > wake_threshold, 1, 'first');

% Use individual cutoffs for lift and drag
C_l_cutoff_idx = C_l_nonzero_idx;
C_d_cutoff_idx = C_d_nonzero_idx;
C_l_cutoff_angle = AoA_total(C_l_cutoff_idx);
C_d_cutoff_angle = AoA_total(C_d_cutoff_idx);

fprintf('C_l coefficients become negligible after %.1f degrees (index %d)\n', C_l_cutoff_angle, C_l_cutoff_idx);
fprintf('C_d coefficients become negligible after %.1f degrees (index %d)\n', C_d_cutoff_angle, C_d_cutoff_idx);

%% Fit functions for RAM coefficients (full range)
fprintf('Fitting RAM coefficients (full AoA range):\n');
ram_degrees = [3, 5, 7,9];
fits = struct();

% Fit C_l
[fits.C_l, fits.C_l_degree] = fit_best_polynomial(AoA_total(C_l_cutoff_idx:end), C_l_total(C_l_cutoff_idx:end), ram_degrees, 'C_l');
% Fit C_d
[fits.C_d, fits.C_d_degree] = fit_best_polynomial(AoA_total(C_d_cutoff_idx:end), C_d_total(C_d_cutoff_idx:end), ram_degrees, 'C_d');
%% Create comprehensive fitted functions that handle full AoA range
AoA_eval = linspace(-90, 90, 1000);

% coefficients (piecewise: fitted function up to cutoff, then zero)
C_l_fitted_poly = zeros(size(AoA_eval));
C_d_fitted_poly = zeros(size(AoA_eval));
C_l_mask = AoA_eval >= C_l_cutoff_angle;
C_d_mask = AoA_eval >= C_d_cutoff_angle;
C_l_fitted_poly(C_l_mask) = fits.C_l(AoA_eval(C_l_mask));
C_d_fitted_poly(C_d_mask) = fits.C_d(AoA_eval(C_d_mask));

%% Plotting
figure();

% Plot 1: Lift coefficients
subplot(2, 2, 1);
plot(AoA_total, C_l_total, 'bo', 'MarkerSize', 4, 'DisplayName', 'C_l (data)');
hold on;
plot(AoA_eval, C_l_fitted_poly, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C_l fit (poly %d)', fits.C_l_degree));

% Mark cutoff
xline(C_l_cutoff_angle, 'k--', 'Alpha', 0.5, 'DisplayName', 'C_l cutoff');
xlabel('Angle of Attack [deg]');
ylabel('Lift Coefficient [-]');
title('Lift Coefficients');
legend('Location', 'best');
grid on;

% Plot 2: Drag coefficients
subplot(2, 2, 2);
plot(AoA_total, C_d_total, 'bo', 'MarkerSize', 4, 'DisplayName', 'C_d (data)');
hold on;
plot(AoA_eval, C_d_fitted_poly, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('C_d fit (poly %d)', fits.C_d_degree));

% Mark cutoff
xline(C_d_cutoff_angle, 'k--', 'Alpha', 0.5, 'DisplayName', 'C_d cutoff');
xlabel('Angle of Attack [deg]');
ylabel('Drag Coefficient [-]');
title('Drag Coefficients');
legend('Location', 'best');
grid on;

% Plot 3: resiudals for lift fits
subplot(2, 2, 3);
AoA_C_l_mask = AoA_total >= C_l_cutoff_angle;
C_l_pred = fits.C_l(AoA_total(AoA_C_l_mask));
C_l_nz = C_l_total(AoA_C_l_mask);
plot(AoA_total(AoA_C_l_mask), C_l_nz - C_l_pred, 'bo-', 'MarkerSize', 3, 'DisplayName', 'C_l residuals');
xlabel('Angle of Attack [deg]');
ylabel('Residuals [-]');
title('C_l Fitting Residuals');
legend('Location', 'best');
grid on;
yline(0, 'k--', 'Alpha', 0.5);

% Plot 4: Residuals for drag fits
subplot(2, 2, 4);
AoA_C_d_mask = AoA_total >= C_d_cutoff_angle;
C_d_pred = fits.C_d(AoA_total(AoA_C_d_mask));
C_d_nz = C_d_total(AoA_C_d_mask);
plot(AoA_total(AoA_C_d_mask), C_d_nz - C_d_pred, 'ro-', 'MarkerSize', 3, 'DisplayName', 'C_d residuals');

xlabel('Angle of Attack [deg]');
ylabel('Residuals [-]');
title('C_d Fitting Residuals');
legend('Location', 'best');
grid on;
yline(0, 'k--', 'Alpha', 0.5);

sgtitle('Aerodynamic Coefficient Function Fitting Analysis');

%% Generate high-resolution AoA vector
AoA_resolution = 1; % 1 degree resolution
AoA_min = -90;
AoA_max = 90;
AoA_highres = (AoA_min:AoA_resolution:AoA_max)';

fprintf('Generating CSV with AoA resolution of %.1e degrees\n', AoA_resolution);
fprintf('AoA range: %.1f to %.1f degrees\n', AoA_min, AoA_max);
fprintf('Total number of points: %d\n', length(AoA_highres));

%% Evaluate fitted functions at high resolution
fprintf('Evaluating fitted functions...\n');

% Initialize coefficient arrays
C_l_highres = zeros(size(AoA_highres));
C_d_highres = zeros(size(AoA_highres));

% For angles below cutoffs: use original data (interpolated)
% For angles above cutoffs: use fitted functions

% Create masks for different regions
below_C_l_cutoff = AoA_highres < C_l_cutoff_angle;
below_C_d_cutoff = AoA_highres < C_d_cutoff_angle;
above_C_l_cutoff = AoA_highres >= C_l_cutoff_angle;
above_C_d_cutoff = AoA_highres >= C_d_cutoff_angle;

% For C_l: use original data below cutoff, fitted function above
if any(below_C_l_cutoff)
    % Interpolate original data for angles below cutoff
    original_mask = AoA_total < C_l_cutoff_angle;
    C_l_highres(below_C_l_cutoff) = interp1(AoA_total(original_mask), C_l_total(original_mask), ...
                                           AoA_highres(below_C_l_cutoff), 'linear', 'extrap');
end
if any(above_C_l_cutoff)
    % Use fitted function for angles above cutoff
    C_l_highres(above_C_l_cutoff) = fits.C_l(AoA_highres(above_C_l_cutoff));
end

% For C_d: use original data below cutoff, fitted function above
if any(below_C_d_cutoff)
    % Interpolate original data for angles below cutoff
    original_mask = AoA_total < C_d_cutoff_angle;
    C_d_highres(below_C_d_cutoff) = interp1(AoA_total(original_mask), C_d_total(original_mask), ...
                                           AoA_highres(below_C_d_cutoff), 'linear', 'extrap');
end
if any(above_C_d_cutoff)
    % Use fitted function for angles above cutoff
    C_d_highres(above_C_d_cutoff) = fits.C_d(AoA_highres(above_C_d_cutoff));
end

%make sure 90 AoA has no Lift
C_l_highres(end) = 0;
%% Create output table
fprintf('Creating output table...\n');

output_table = table(AoA_highres, C_l_highres, C_d_highres, ...
                     'VariableNames', {'AoA', 'C_l', 'C_d'});

%% Write to CSV file

fprintf('Writing to file: %s\n', output_filename);

writetable(output_table, output_filename);
%% Helper function for polynomial fitting
function [best_fit_func, best_degree] = fit_best_polynomial(x, y, degrees, coeff_name)
    best_rmse = inf;
    best_fit_func = [];
    best_degree = degrees(1);
    
    for deg = degrees
        p = polyfit(x, y, deg);
        y_pred = polyval(p, x);
        rmse = sqrt(mean((y - y_pred).^2));
        
        fprintf('%s degree %d: RMSE = %.2e\n', coeff_name, deg, rmse);
        
        if rmse < best_rmse
            best_rmse = rmse;
            best_fit_func = @(x_eval) polyval(p, x_eval);
            best_degree = deg;
        end
    end
    
    fprintf('Best fit for %s: degree %d (RMSE = %.2e)\n', coeff_name, best_degree, best_rmse);
end
