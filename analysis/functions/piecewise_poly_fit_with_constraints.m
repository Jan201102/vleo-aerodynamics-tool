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