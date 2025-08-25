function damping = calculate_aerodynamic_damping(method, delta__rad_per_s, axis_dir_B, component, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions)
%% calculate_aerodynamic_damping - calculate the aerodynamic damping using numerical differentiation
%
% damping = calculate_aerodynamic_damping(method,...
%                                             delta__rad_per_s,...
%                                             axis_dir_B,...
%                                             component,...
%                                             bodies,...
%                                             bodies_rotation_angles__rad,...
%                                             model,...
%                                             lut_data,...
%                                             environment_definitions)
% This function calculates the aerodynamic damping for a given method of numerical differentiation.
%
% Inputs:
%   method: string specifying the numerical differentiation method to use
%           Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward' 
%   delta__rad_per_s: step size for numerical differentiation
%   axis_dir_B: 3x1 vector specifying the rotation axis direction in body frame
%   component: index of the torque component to differentiate (1 = x, 2 = y, 3 = z)
%   bodies: cell array of body structures containing aerodynamic properties
%   bodies_rotation_angles__rad: 1xN vector of rotation angles for each body in radians
%   model: integer specifying the aerodynamic model to use
%   lut_data: gridded interpolant for LUT data
%   environment_definitions: struct containing environment parameters
%
% Outputs:
%   damping: scalar value representing the aerodynamic damping for the specified method
%   calculated at the specified rotation axis and component
    arguments
        method (1,1) string {mustBeMember(method, {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'})}
        delta__rad_per_s (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        axis_dir_B (3,1) {mustBeNumeric, mustBeReal}; % Add rotation axis parameter
        component (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBeGreaterThanOrEqual(component, 1), mustBeLessThanOrEqual(component, 3)};
        bodies (1,:) cell {mustBeNonempty};  % Cell array of body structures
        bodies_rotation_angles__rad (1,:) {mustBeNumeric, mustBeReal};  % Rotation angles for each body in radians
        model (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive};  % Aerodynamic model index
        lut_data ; % lookup table data for aerodynamic coefficients
        environment_definitions struct; % Environment parameters
    end
    
    % Create function handle with parameterized rotation axis
    torque_func = @(rotational_velocity_magnitude) compute_torque_at_rotational_velocity(rotational_velocity_magnitude, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions);
    
    % Use generic numerical differentiation function (derivative at rotational velocity = 0)
    damping = numerical_differentiation(method, delta__rad_per_s, torque_func, 0, component);
end

function torque = compute_torque_at_rotational_velocity(rotational_velocity_magnitude, axis_direction, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions)
    import vleo_aerodynamics_core.*
    
    attitude_quaternion_BI = [1; 0; 0; 0]; % Identity quaternion
    
    % Apply rotational velocity around the specified axis
    rotational_velocity_vector = rotational_velocity_magnitude * axis_direction;
    
    [~, torque, ~, ~] = vleoAerodynamics(...
        attitude_quaternion_BI,...
        rotational_velocity_vector,...
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
end