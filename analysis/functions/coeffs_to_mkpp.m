
function pp = coeffs_to_mkpp(coeffs, knots, degrees)
%COEFFS_TO_MKPP Convert piecewise polynomial coefficients to MATLAB pp structure
%
%   pp = COEFFS_TO_MKPP(coeffs, knots, degrees)
%
%   Converts global coordinate piecewise polynomial coefficients from 
%   piecewise_poly_fit_with_constraints to a MATLAB piecewise polynomial
%   structure compatible with ppval.
%
% Inputs:
%   coeffs  - Coefficient vector from piecewise_poly_fit_with_constraints
%   knots   - Vector of knot points (same as used in fitting)
%   degrees - Vector of polynomial degrees for each segment
%
% Outputs:
%   pp - MATLAB piecewise polynomial structure compatible with ppval

    %% Initialize Parameters
    n_segments = length(knots) + 1;
    coeff_counts = degrees + 1;
    
    % Create break points (including endpoints)
    % Note: Endpoints may need adjustment based on your data range
    breaks = [-90, knots, 90];
    
    % Initialize coefficient matrix for mkpp format
    max_degree = max(degrees);
    coeff_matrix = zeros(n_segments, max_degree + 1);
    
    %% Convert Each Segment from Global to Local Coordinates
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