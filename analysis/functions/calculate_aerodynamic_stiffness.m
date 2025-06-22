function stiffness = calculate_aerodynamic_stiffness(method, delta__rad, axis_dir_B, component, bodies, bodies_rotation_angles__rad, model, lut_data)
%% calculate_aerodynamic_stiffness - calculate the aerodynamic stiffness using numerical differentiation
%
% stiffness = calculate_aerodynamic_stiffness(method,...
%                                             delta,...
%                                             axis_dir,...
%                                             component,...
%                                             bodies,...
%                                             bodies_rotation_angles__rad,...
%                                             model,...
%                                             lut_data)
% This function calculates the aerodynamic stiffness for a given method of numerical differentiation.
%
% Inputs:
%   method: string specifying the numerical differentiation method to use
%           Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward' 
%   delta: step size for numerical differentiation
%   axis_dir_B: 3x1 vector specifying the rotation axis in the body frame
%   component: index of the torque component to differentiate (1 = x, 2 = y, 3 = z)
%   bodies: cell array of body structures containing aerodynamic properties
%   bodies_rotation_angles__rad: 1xN vector of rotation angles for each body in radians
%   model: integer specifying the aerodynamic model to use
%   lut_data: lookup table data for aerodynamic coefficients
%
% Outputs:
%   stiffness: scalar value representing the aerodynamic stiffness for the specified method
%   calculated at the specified rotation axis and component
    arguments
        method (1,1) string {mustBeMember(method, {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'})}
        delta__rad (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        axis_dir_B (3,1) {mustBeNumeric, mustBeReal};
        component (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBeGreaterThanOrEqual(component, 1), mustBeLessThanOrEqual(component, 3)};
        bodies (1,:) cell {mustBeNonempty};  % Cell array of body structures
        bodies_rotation_angles__rad (1,:) {mustBeNumeric, mustBeReal};  % Rotation angles for each body in radians
        model (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive};  % Aerodynamic model index
        lut_data (:,5) {mustBeNumeric, mustBeReal}; % lookup table data for aerodynamic coefficients
    end

    switch method
        case 'central'
            % Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
            torque_plus = compute_aerodynamics_at_angle(delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_minus = compute_aerodynamics_at_angle(-delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            stiffness = (torque_plus(component) - torque_minus(component)) / (2 * delta__rad);
            
        case 'five_point_stencil'
            % Five-point stencil: f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
            torque_1 = compute_aerodynamics_at_angle(-2*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_2 = compute_aerodynamics_at_angle(-delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_3 = compute_aerodynamics_at_angle(delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_4 = compute_aerodynamics_at_angle(2*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            stiffness = (-torque_4(component) + 8*torque_3(component) - 8*torque_2(component) + torque_1(component)) / (12 * delta__rad);
            
        case 'seven_point_stencil'
            % Seven-point stencil: f'(x) ≈ (-f(x-3h) + 9f(x-2h) - 45f(x-h) + 45f(x+h)  9f(x+2h) + f(x+3h)) / (60h)
            torque_1 = compute_aerodynamics_at_angle(-3*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_2 = compute_aerodynamics_at_angle(-2*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_3 = compute_aerodynamics_at_angle(-delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_4 = compute_aerodynamics_at_angle(delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_5 = compute_aerodynamics_at_angle(2*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_6 = compute_aerodynamics_at_angle(3*delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            stiffness = (-torque_1(component) + 9*torque_2(component) - 45*torque_3(component) + 45*torque_4(component) - 9*torque_5(component) + torque_6(component)) / (60 * delta__rad);
            
        case 'forward'
            % Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
            torque_0 = compute_aerodynamics_at_angle(0, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_plus = compute_aerodynamics_at_angle(delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            stiffness = (torque_plus(component) - torque_0(component)) / delta__rad;
            
        case 'backward'
            % Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
            torque_0 = compute_aerodynamics_at_angle(0, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_minus = compute_aerodynamics_at_angle(-delta__rad, axis_dir_B, bodies, bodies_rotation_angles__rad, model, lut_data);
            stiffness = (torque_0(component) - torque_minus(component)) / delta__rad;
            
        otherwise
            error('Unknown derivation method: %s', method);
    end
end
%% Helper Functions

function torque = compute_aerodynamics_at_angle(attitude_angle__rad, axis_direction, bodies, bodies_rotation_angles__rad, model, lut_data)
    import vleo_aerodynamics_core.*
    % Compute aerodynamic torque at a specific attitude angle
    attitude_quarternion_BI = [cos(attitude_angle__rad/2); sin(attitude_angle__rad/2) * axis_direction];
    
    % Get environment variables from caller workspace
    rotational_velocity_BI_B__rad_per_s = evalin('base', 'rotational_velocity_BI_B__rad_per_s');
    velocity_I_I__m_per_s = evalin('base', 'velocity_I_I__m_per_s');
    wind_velocity_I_I__m_per_s = evalin('base', 'wind_velocity_I_I__m_per_s');
    density__kg_per_m3 = evalin('base', 'density__kg_per_m3');
    temperature__K = evalin('base', 'temperature__K');
    particles_mass__kg = evalin('base', 'particles_mass__kg');
    temperature_ratio_method = evalin('base', 'temperature_ratio_method');
    
    [~, torque,~,~] = vleoAerodynamics(...
        attitude_quarternion_BI,...
        rotational_velocity_BI_B__rad_per_s,...
        velocity_I_I__m_per_s,...
        wind_velocity_I_I__m_per_s,...
        density__kg_per_m3,...
        temperature__K,...
        particles_mass__kg,...
        bodies,...
        bodies_rotation_angles__rad,...
        temperature_ratio_method,...
        model,...
        lut_data);
end