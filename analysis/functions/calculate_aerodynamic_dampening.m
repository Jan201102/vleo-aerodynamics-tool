function dampening = calculate_aerodynamic_dampening(method, delta__rad_per_s, component, bodies, bodies_rotation_angles__rad,model, lut_data)
%% calculate_aerodynamic_dampening - calculate the aerodynamic dampening using numerical differentiation
%
% stiffness = calculate_aerodynamic_stiffness(method,...
%                                             delta__rad_per_s,...
%                                             axis_dir,...
%                                             component,...
%                                             bodies,...
%                                             bodies_rotation_angles__rad,...
%                                             model,...
%                                             lut_data)
% This function calculates the aerodynamic dampening for a given method of numerical differentiation.
%
% Inputs:
%   method: string specifying the numerical differentiation method to use
%           Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward' 
%   delta__rad_per_s: step size for numerical differentiation
%   component: index of the torque component to differentiate (1 = x, 2 = y, 3 = z)
%   bodies: cell array of body structures containing aerodynamic properties
%   bodies_rotation_angles__rad: 1xN vector of rotation angles for each body in radians
%   model: integer specifying the aerodynamic model to use
%   lut_data: lookup table data for aerodynamic coefficients
%
% Outputs:
%   dampening: scalar value representing the aerodynamic dampening for the specified method
%   calculated at the specified rotation axis and component
    arguments
        method (1,1) string {mustBeMember(method, {'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'})}
        delta__rad_per_s (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        component (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBeGreaterThanOrEqual(component, 1), mustBeLessThanOrEqual(component, 3)};
        bodies (1,:) cell {mustBeNonempty};  % Cell array of body structures
        bodies_rotation_angles__rad (1,:) {mustBeNumeric, mustBeReal};  % Rotation angles for each body in radians
        model (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive};  % Aerodynamic model index
        lut_data (:,5) {mustBeNumeric, mustBeReal}; % lookup table data for aerodynamic coefficients
    end
    switch method
        case 'central'
            % Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
            torque_plus = compute_aerodynamics_at_rotational_velocity(delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_minus = compute_aerodynamics_at_rotational_velocity(-delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            dampening = (torque_plus(component) - torque_minus(component)) / (2 * delta__rad_per_s);
            
        case 'five_point_stencil'
            % Five-point stencil: f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
            torque_1 = compute_aerodynamics_at_rotational_velocity(-2*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_2 = compute_aerodynamics_at_rotational_velocity(-delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_3 = compute_aerodynamics_at_rotational_velocity(delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_4 = compute_aerodynamics_at_rotational_velocity(2*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            dampening = (-torque_4(component) + 8*torque_3(component) - 8*torque_2(component) + torque_1(component)) / (12 * delta__rad_per_s);
            
        case 'seven_point_stencil'
            % Seven-point stencil: f'(x) ≈ (-f(x-3h) + 9f(x-2h) - 45f(x-h) + 45f(x+h) + 9f(x+2h) + f(x+3h)) / (60h)
            torque_1 = compute_aerodynamics_at_rotational_velocity(-3*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_2 = compute_aerodynamics_at_rotational_velocity(-2*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_3 = compute_aerodynamics_at_rotational_velocity(-delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_4 = compute_aerodynamics_at_rotational_velocity(delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_5 = compute_aerodynamics_at_rotational_velocity(2*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_6 = compute_aerodynamics_at_rotational_velocity(3*delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            dampening = (-torque_6(component) + 9*torque_5(component) - 45*torque_4(component) + 45*torque_3(component) + 9*torque_2(component) + torque_1(component)) / (60 * delta__rad_per_s);
        case 'forward'
            % Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
            torque_plus = compute_aerodynamics_at_rotational_velocity(delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_current = compute_aerodynamics_at_rotational_velocity(0, bodies, bodies_rotation_angles__rad, model, lut_data);
            dampening = (torque_plus(component) - torque_current(component)) / delta__rad_per_s;
        case 'backward'
            % Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
            torque_current = compute_aerodynamics_at_rotational_velocity(0, bodies, bodies_rotation_angles__rad, model, lut_data);
            torque_minus = compute_aerodynamics_at_rotational_velocity(-delta__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data);
            dampening = (torque_current(component) - torque_minus(component)) / delta__rad_per_s;
    end
end

%% Helper Functions
function torque = compute_aerodynamics_at_rotational_velocity(rotational_velocity_BI_B__rad_per_s, bodies, bodies_rotation_angles__rad, model, lut_data)
    import vleo_aerodynamics_core.*
    % Compute aerodynamic torque at a specific attitude angle
    attitude_quarternion_BI = [1; 0; 0; 0]; % No rotation, identity quaternion
    
    % Get environment variables from caller workspace
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