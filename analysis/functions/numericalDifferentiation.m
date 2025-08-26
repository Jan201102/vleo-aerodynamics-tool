function derivative = numerical_differentiation(method, delta, func_handle, x0, component)
%% numerical_differentiation - Generic numerical differentiation function
%
% derivative = numerical_differentiation(method, delta, func_handle, x0, component)
%
% This function calculates the numerical derivative of a function using
% various finite difference methods.
%
% Inputs:
%   method: string specifying the numerical differentiation method to use
%           Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward' 
%   delta: step size for numerical differentiation
%   func_handle: function handle to differentiate
%   x0: point at which to evaluate the derivative
%   component: (optional) component index for vector-valued functions
%
% Outputs:
%   derivative: scalar or vector value representing the numerical derivative
    arguments
        method (1,1) string {mustBeMember(method, {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'})}
        delta (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        func_handle function_handle;
        x0 (1,1) {mustBeNumeric, mustBeReal};
        component (1,1) {mustBeNumeric, mustBeReal, mustBeInteger} = [];
    end
    
    % Helper function to extract component if needed
    extract_component = @(vec) if_component_specified(vec, component);
    
    switch method
        case 'central'
            % Central difference: f'(x) ~= (f(x+h) - f(x-h)) / (2h)
            f_plus = extract_component(func_handle(x0 + delta));
            f_minus = extract_component(func_handle(x0 - delta));
            derivative = (f_plus - f_minus) / (2 * delta);
            
        case 'five_point_stencil'
            % Five-point stencil: f'(x) ~= (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
            f_1 = extract_component(func_handle(x0 - 2*delta));
            f_2 = extract_component(func_handle(x0 - delta));
            f_3 = extract_component(func_handle(x0 + delta));
            f_4 = extract_component(func_handle(x0 + 2*delta));
            derivative = (-f_4 + 8*f_3 - 8*f_2 + f_1) / (12 * delta);
            
        case 'seven_point_stencil'
            % Seven-point stencil: f'(x) ~= (-f(x-3h) + 9f(x-2h) - 45f(x-h) + 45f(x+h) - 9f(x+2h) + f(x+3h)) / (60h)
            f_1 = extract_component(func_handle(x0 - 3*delta));
            f_2 = extract_component(func_handle(x0 - 2*delta));
            f_3 = extract_component(func_handle(x0 - delta));
            f_4 = extract_component(func_handle(x0 + delta));
            f_5 = extract_component(func_handle(x0 + 2*delta));
            f_6 = extract_component(func_handle(x0 + 3*delta));
            derivative = (-f_1 + 9*f_2 - 45*f_3 + 45*f_4 - 9*f_5 + f_6) / (60 * delta);
            
        case 'forward'
            % Forward difference: f'(x) ~= (f(x+h) - f(x)) / h
            f_0 = extract_component(func_handle(x0));
            f_plus = extract_component(func_handle(x0 + delta));
            derivative = (f_plus - f_0) / delta;
            
        case 'backward'
            % Backward difference: f'(x) ~= (f(x) - f(x-h)) / h
            f_0 = extract_component(func_handle(x0));
            f_minus = extract_component(func_handle(x0 - delta));
            derivative = (f_0 - f_minus) / delta;
            
        otherwise
            error('Unknown differentiation method: %s', method);
    end
end

function result = if_component_specified(vec, component)
    if isempty(component)
        result = vec;
    else
        result = vec(component);
    end
end
