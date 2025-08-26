function stiffness = calculate_aerodynamic_stiffness(method, alpha_zero, delta__rad, axis_dir_B, component, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions)
%% calculate_aerodynamic_stiffness - calculate the aerodynamic stiffness using numerical differentiation
%
% stiffness = calculate_aerodynamic_stiffness(method,...
%                                             alpha_zero,...
%                                             delta__rad,...
%                                             axis_dir_B,...
%                                             component,...
%                                             bodies,...
%                                             bodies_rotation_angles__rad,...
%                                             model,...
%                                             lut_data,...
%                                             environment_definitions)
% This function calculates the aerodynamic stiffness for a given method of numerical differentiation.
%
% Inputs:
%   method: string specifying the numerical differentiation method to use
%           Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward' 
%   alpha_zero: angle the stiffness is calculated at, in radians
%   delta__rad: step size for numerical differentiation
%   axis_dir_B: 3x1 vector specifying the rotation axis in the body frame
%   component: index of the torque component to differentiate (1 = x, 2 = y, 3 = z)
%   bodies: cell array of body structures containing aerodynamic properties
%   bodies_rotation_angles__rad: 1xN vector of rotation angles for each body in radians
%   model: integer specifying the aerodynamic model to use
%   lut_data: lookup table data of type 'griddedInterpolant' for aerodynamic coefficients
%   environment_definitions: struct containing environment parameters
%
% Outputs:
%   stiffness: scalar value representing the aerodynamic stiffness for the specified method
%   calculated at the specified rotation axis and component
    arguments
        method (1,1) string {mustBeMember(method, {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'})}
        alpha_zero (1,1) {mustBeNumeric, mustBeReal};
        delta__rad (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        axis_dir_B (3,1) {mustBeNumeric, mustBeReal};
        component (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBeGreaterThanOrEqual(component, 1), mustBeLessThanOrEqual(component, 3)};
        bodies (1,:) cell {mustBeNonempty};  % Cell array of body structures
        bodies_rotation_angles__rad (1,:) {mustBeNumeric, mustBeReal};  % Rotation angles for each body in radians
        model (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive};  % Aerodynamic model index
        lut_data; % lookup table data for aerodynamic coefficients
        environment_definitions struct; % Environment parameters
    end

    % Create function handle for torque computation at different angles
    torque_func = @(attitude_angle__rad) compute_torque_at_angle(attitude_angle__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions);
    
    % Use generic numerical differentiation function
    stiffness = numerical_differentiation(method, delta__rad, torque_func, alpha_zero, component);
    
end

function torque = compute_torque_at_angle(attitude_angle__rad, axis_direction, bodies, bodies_rotation_angles__rad, model, lut_data, environment_definitions)
    import vleo_aerodynamics_core.*
    % Compute aerodynamic torque at a specific attitude angle
    attitude_quaternion_BI = [cos(attitude_angle__rad/2); sin(attitude_angle__rad/2) * axis_direction];
    
    [~, torque, ~, ~] = vleoAerodynamics(...
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
end