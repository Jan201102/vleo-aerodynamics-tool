%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
clear;

%load model data
bodies = load_model();

%load environment data
run("environment_definitions.m");

%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'cl_cd_cVAE_A01_flat_and_bird.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% Derivation Configuration
derivation_method = 'seven_point_stencil';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_alpha = pi/1800000;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% show bodies
showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);

%% calculate aerodynamic stiffness
num_angles = 101;
control_surface_angles__rad = linspace(pi/180, pi/2-pi/180, num_angles);
aero_stiffness_pitch = nan(num_angles,2);

for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = ones(1,5) * current_angle;
        
        % Calculate numerical derivative using configurable method
        aero_stiffness_pitch(i,model) = calculate_aerodynamic_derivative(...
            derivation_method, delta_alpha, axis_direction, torque_component, ...
            bodies, bodies_rotation_angles__rad, model, lut);
    end
end

%% plot aerodynamic stiffness
figure;
tl = tiledlayout('flow');
title(tl, 'Aerodynamic Stiffness Calculation');

ax1 = nexttile;
plot(ax1,control_surface_angles__rad,aero_stiffness_pitch(:,1), 'LineWidth', 2);
grid on;
title(' Sentman Model');
xlabel('Control Surface Angle [rad]');
ylabel('Aerodynamic Stiffness [Nm/rad]');

ax2 = nexttile;
plot(ax2,control_surface_angles__rad, aero_stiffness_pitch(:,2), 'LineWidth', 2);
grid on;
title('IRS Model');
xlabel('Control Surface Angle [rad]');
ylabel('Aerodynamic Stiffness [Nm/rad]');

%% Helper Functions

function torque = compute_aerodynamics_at_angle(attitude_angle, axis_direction, bodies, bodies_rotation_angles__rad, model, lut)
    import vleo_aerodynamics_core.*
    % Compute aerodynamic torque at a specific attitude angle
    attitude_quarternion_BI = [cos(attitude_angle/2); sin(attitude_angle/2) * axis_direction];
    
    % Get environment variables from caller workspace
    rotational_velocity_BI_B__rad_per_s = evalin('base', 'rotational_velocity_BI_B__rad_per_s');
    velocity_I_I__m_per_s = evalin('base', 'velocity_I_I__m_per_s');
    wind_velocity_I_I__m_per_s = evalin('base', 'wind_velocity_I_I__m_per_s');
    density__kg_per_m3 = evalin('base', 'density__kg_per_m3');
    temperature__K = evalin('base', 'temperature__K');
    particles_mass__kg = evalin('base', 'particles_mass__kg');
    temperature_ratio_method = evalin('base', 'temperature_ratio_method');
    
    [~, torque] = vleoAerodynamics(...
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
        lut);
end

function derivative = calculate_aerodynamic_derivative(method, delta, axis_dir, component, bodies, bodies_rotation_angles__rad, model, lut)
    % Calculate numerical derivative using specified method
    
    switch method
        case 'central'
            % Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
            torque_plus = compute_aerodynamics_at_angle(delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_minus = compute_aerodynamics_at_angle(-delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            derivative = (torque_plus(component) - torque_minus(component)) / (2 * delta);
            
        case 'five_point_stencil'
            % Five-point stencil: f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
            torque_1 = compute_aerodynamics_at_angle(-2*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_2 = compute_aerodynamics_at_angle(-delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_3 = compute_aerodynamics_at_angle(delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_4 = compute_aerodynamics_at_angle(2*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            derivative = (-torque_4(component) + 8*torque_3(component) - 8*torque_2(component) + torque_1(component)) / (12 * delta);
            
        case 'seven_point_stencil'
            % Seven-point stencil: f'(x) ≈ (-f(x-3h) + 9f(x-2h) - 45f(x-h) + 45f(x+h)  9f(x+2h) + f(x+3h)) / (60h)
            torque_1 = compute_aerodynamics_at_angle(-3*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_2 = compute_aerodynamics_at_angle(-2*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_3 = compute_aerodynamics_at_angle(-delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_4 = compute_aerodynamics_at_angle(delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_5 = compute_aerodynamics_at_angle(2*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_6 = compute_aerodynamics_at_angle(3*delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            derivative = (-torque_1(component) + 9*torque_2(component) - 45*torque_3(component) + 45*torque_4(component) - 9*torque_5(component) + torque_6(component)) / (60 * delta);
            
        case 'forward'
            % Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
            torque_0 = compute_aerodynamics_at_angle(0, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_plus = compute_aerodynamics_at_angle(delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            derivative = (torque_plus(component) - torque_0(component)) / delta;
            
        case 'backward'
            % Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
            torque_0 = compute_aerodynamics_at_angle(0, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            torque_minus = compute_aerodynamics_at_angle(-delta, axis_dir, bodies, bodies_rotation_angles__rad, model, lut);
            derivative = (torque_0(component) - torque_minus(component)) / delta;
            
        otherwise
            error('Unknown derivation method: %s', method);
    end
end