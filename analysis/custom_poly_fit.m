%% custom poly fit
clear; clc; close all;
csv_file_in = 'aerodynamic_coefficients_panel_method_sentman.csv';
data = readtable(csv_file_in);
% Extract data
AoA = data.AoA;
C_l_ram = data.C_l_ram;
%ensure 90 lift value is 0
C_l_ram(end) = 0;
C_d_ram = data.C_d_ram;
C_l_wake = data.C_l_wake;
C_d_wake = data.C_d_wake;
% and reorder concatenate data
AoA_total = [flip(-AoA(2:end));AoA];
C_l_total = [flip(-C_l_wake(2:end));C_l_ram];
C_d_total = [flip(C_d_wake(2:end));C_d_ram];
data = [AoA_total C_l_total C_d_total];

%% Define fitting configuration
% Boundaries define the ranges: [range_min, boundary_1_2, boundary_2_3, range_max]
boundaries = [-24, 0];

% Polynomial degrees for each range
poly_degrees = [0, 9, 9];
poly_coeffs = poly_degrees + 1; % Number of coefficients for each polynomial
% Constraint points for C_l (cell array)
constraint_points_Cl = [90,0] ;          % Range 3 constraints
global_piecewise_coeff_cl = piecewise_poly_fit_with_constraints(AoA_total, C_l_total, boundaries, poly_degrees, constraint_points_Cl,1);
global_piecewise_coeff_cd = piecewise_poly_fit_with_constraints(AoA_total, C_d_total, boundaries, poly_degrees, [],1);
pp_Cl = coeffs_to_mkpp(global_piecewise_coeff_cl, boundaries, poly_degrees);
pp_Cd = coeffs_to_mkpp(global_piecewise_coeff_cd, boundaries, poly_degrees);
save("aerodynamic_coefficients_panel_method_sentman_poly.mat","pp_Cd","pp_Cl");
%% Test and plot the results
AoA_test = linspace(-90, 90, 1000);
Cl_combined = ppval(pp_Cl, AoA_test);
Cd_combined = ppval(pp_Cd, AoA_test);

% Plot C_l
figure;
subplot(2,1,1);
hold on;
scatter(AoA_total, C_l_total, 'ko', 'DisplayName', 'Original C_l Data');
plot(AoA_test, Cl_combined, 'm-', 'LineWidth', 3, 'DisplayName', 'Combined Piecewise Polynomial');
xlabel('Angle of Attack (degrees)');
ylabel('C_l');
title('C_l Piecewise Polynomial Fit');
legend('Location', 'best');
grid on;

% Plot C_d
subplot(2,1,2);
hold on;
scatter(AoA_total, C_d_total, 'ko', 'DisplayName', 'Original C_d Data');
plot(AoA_test, Cd_combined, 'm-', 'LineWidth', 3, 'DisplayName', 'Combined Piecewise Polynomial');

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

%% piecewise polynomial fitting with constraints
function coeffs = piecewise_poly_fit_with_constraints(x, y, knots, degrees,constraints, continuity_degree)
    % Fit a piecewise polynomial with constraints
    % Inputs:
    %   x: vector of x data points
    %   y: vector of y data points
    %   knots 1xN: vector of knot points defining the piecewise segments
    %   degrees 1x(N+1): vector of polynomial degrees for each segment
    %   constraints: Mx2 matrix where each row is [constraint_x, constraint_y]
    %   continuity_degree: degree of continuity (0: C0 continuity, 1: C1 continuity, 2: C2 continuity)
    % Outputs:
    %   coeffs: coefficients of the piecewise polynomial

    % Ensure knots are sorted
    knots = sort(knots);
    n_segments = length(knots) + 1;
    if length(degrees) ~= n_segments
        error('Degrees must match the number of segments defined by knots.');
    end
    
    % Initialize the coefficient matrix and right-hand side vector
    n_constraints = size(constraints, 1);
    coeff_counts = degrees+1;
    n_coeffs = sum(coeff_counts);
    A = zeros(length(y), n_coeffs);
    
    % Create the Vandermonde matrix for each segment
    for i = 1:n_segments
        if i == 1
            segment_idx = x <= knots(i);
        elseif i == n_segments
            segment_idx = x > knots(i-1);
        else
            segment_idx = x > knots(i-1) & x <= knots(i);
        end
        
        % Create Vandermonde matrix for the current segment
        for k = 1:coeff_counts(i)
            col_idx = sum(coeff_counts(1:i-1)) + k;
            A(segment_idx, col_idx) = x(segment_idx).^(k-1);
        end
    end

    %create constraint matrix
    C = zeros(n_constraints+(continuity_degree+1)*length(knots), n_coeffs);
    d = zeros(n_constraints+(continuity_degree+1)*length(knots), 1);
    
    % Fill the constraint matrix
    for i = 1:n_constraints
        x_const = constraints(i, 1);
        %determine which segment the constraint belongs to
        segment_idx = find(knots >= x_const, 1);
        if isempty(segment_idx)
            segment_idx = n_segments; % last segment if x_const is beyond last knot
        end
        % Fill the constraint matrix
        for k = 1:coeff_counts(segment_idx)
            col_idx = sum(coeff_counts(1:segment_idx-1)) + k;
            C(i, col_idx) = x_const^(k-1);
        end
        d(i) = constraints(i, 2);
    end

    % Add continuity constraints
    for i = 1:continuity_degree+1
        for k = 1:length(knots)
            %left side of the knot
            for c = 1:coeff_counts(k)
                if c >=i
                    col_idx_left = sum(coeff_counts(1:k-1)) + c;
                    C(n_constraints + (i-1)*length(knots) + k, col_idx_left) = factorial(c-1)/factorial(c-i)*knots(k)^(c-i);
                end
            end
            %right side of the knot
            for c = 1:coeff_counts(k+1)
                if c >=i
                    col_idx_right = sum(coeff_counts(1:k)) + c;
                    C(n_constraints + (i-1)*length(knots) + k, col_idx_right) = -factorial(c-1)/factorial(c-i)*knots(k)^(c-i);
                end
            end
        end
    end

    % Solve the system of equations
    % [A'*A  C'] [p]   [A'*y]
    % [C     0 ] [λ] = [d   ]
    ATA = A' * A;
    ATy = A' * y;
    % Build the augmented system
    augmented_matrix = [ATA, C'; C, zeros(size(C, 1))];
    rhs = [ATy; d];
    % Solve the system
    solution = augmented_matrix \ rhs;
    % Extract polynomial coefficients
    coeffs = solution(1:n_coeffs);

end
function pp = coeffs_to_mkpp(coeffs, knots, degrees)
    % Convert global coordinate piecewise polynomial coefficients to mkpp object
    % Inputs:
    %   coeffs: coefficient vector from piecewise_poly_fit_with_constraints
    %   knots: vector of knot points (same as used in fitting)
    %   degrees: vector of polynomial degrees for each segment
    % Output:
    %   pp: MATLAB piecewise polynomial structure compatible with ppval
    
    % Number of segments
    n_segments = length(knots) + 1;
    coeff_counts = degrees + 1;
    
    % Create break points (including endpoints)
    % Assume the first segment starts at the minimum x value in your data
    % and the last segment extends beyond the last knot
    breaks = [-90, knots, 90];  % You may need to adjust the endpoints
    
    % Initialize coefficient matrix for mkpp
    max_degree = max(degrees);
    coeff_matrix = zeros(n_segments, max_degree + 1);
    
    % Convert each segment from global to local coordinates
    for i = 1:n_segments
        % Extract coefficients for this segment (in ascending power order)
        start_idx = sum(coeff_counts(1:i-1)) + 1;
        end_idx = sum(coeff_counts(1:i));
        segment_coeffs = coeffs(start_idx:end_idx);
        
        % Get the left endpoint of this segment for coordinate transformation
        if i == 1
            x_left = breaks(i+1);  % Use first knot as reference
        else
            x_left = breaks(i);    % Use the knot at the left boundary
        end
        
        % Convert from global coordinates p(x) to local coordinates p(x-x_left)
        % If p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
        % Then p(x_left + t) = sum over k of (ak * sum over j of C(k,j) * x_left^(k-j) * t^j)
        
        local_coeffs = zeros(1, degrees(i) + 1);
        
        % Transform using binomial expansion
        for k = 1:length(segment_coeffs)  % k is the power in global coordinates (0-indexed as k-1)
            global_power = k - 1;
            global_coeff = segment_coeffs(k);
            
            % Expand (x_left + t)^global_power using binomial theorem
            for j = 0:global_power
                local_power = j;
                binomial_coeff = nchoosek(global_power, j);
                x_left_power = global_power - j;
                
                if local_power + 1 <= length(local_coeffs)
                    local_coeffs(local_power + 1) = local_coeffs(local_power + 1) + ...
                        global_coeff * binomial_coeff * (x_left^x_left_power);
                end
            end
        end
        
        % Store coefficients in descending power order for mkpp (highest power first)
        local_coeffs_desc = flip(local_coeffs);
        
        % Pad with zeros if necessary to match max_degree
        if length(local_coeffs_desc) < max_degree + 1
            padded_coeffs = [zeros(1, max_degree + 1 - length(local_coeffs_desc)), local_coeffs_desc];
        else
            padded_coeffs = local_coeffs_desc;
        end
        
        coeff_matrix(i, :) = padded_coeffs;
    end
    
    % Create the piecewise polynomial structure
    pp = mkpp(breaks, coeff_matrix);
end