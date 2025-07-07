%% AERODYNAMIC STIFFNESS
% calculate the aerodynamic stiffness for a satellite for different control surface configurations
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLegendFontName', 'Times New Roman');
%load environment data
run("environment_definitions.m");

%% load lut data
lut_data = load_lut("aerodynamic_coefficients_panel_method.csv");

bodies = load_from_gmsh(energy_accommodation,surface_temperature__K);
showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);
num_bodies = 5;
rotation_face_index = [2,3,4,5];
x_label = "control surface angle [°]";
I_yy = 0.0075833; % Moment of inertia around y-axis

%% Derivation Configuration
derivation_method = 'central';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_q_rad_per_s = 1e2;                  % Step size for numerical differentiation
delta_alpha = 1e-4;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% calculate aerodynamic stiffness and dampening
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
ol_parameters = zeros(num_angles,2,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        C_M_da = calculate_aerodynamic_dampening(derivation_method, ...
                                                delta_q_rad_per_s, ...
                                                torque_component, ...
                                                bodies, ...
                                                bodies_rotation_angles__rad, ...
                                                model, ...
                                                lut_data);
        C_M_a = calculate_aerodynamic_stiffness(derivation_method, ...
                                                delta_alpha, ...
                                                axis_direction,...
                                                torque_component, ...
                                                bodies, ...
                                                bodies_rotation_angles__rad, ...
                                                model, ...
                                                lut_data);
        [omega_0,zeta] = get_omega_zeta(I_yy,C_M_a,C_M_da);
        ol_parameters(i,model,1) = omega_0;
        ol_parameters(i,model,2) = zeta;
        fprintf('calcluated point %d of %d\n',i+(model-1)*num_angles,2*num_angles);
    end
end
%% save calculted ol parameters
save("ol_parameters.mat","ol_parameters","control_surface_angles__rad")
%% plot omega_0 and zeta in tiled layout
figure;
t = tiledlayout(2,1);
nexttile;
plot(control_surface_angles__rad * 180 / pi, squeeze(ol_parameters(:,1,1)),"b-");
hold on;
plot(control_surface_angles__rad * 180 / pi, squeeze(ol_parameters(:,2,1)),"r-");
grid on;
xlabel(x_label);
ylabel('Eigenfrequenz \omega_0 [rad/s]');
title('Eigenfrequenz \omega_0 für verschiedene Modelle');
legend("Sentman","IRS","Location","southeast");
nexttile;
plot(control_surface_angles__rad * 180 / pi, squeeze(ol_parameters(:,1,2)),"b-");
hold on;
plot(control_surface_angles__rad * 180 / pi, squeeze(ol_parameters(:,2,2)),"r-");
grid on;
xlabel(x_label);
ylabel('Dämpfungsgrad \zeta [-]');
title('Dämpfungsgrad \zeta für verschiedene Modelle');
legend("Sentman","IRS","Location","southeast");
%%
matlab2tikz('ol_parameters_shuttlecock.tex');
%%
plot_solution_at_angle(deg2rad(45));
%%
function [omega_0,zeta] = get_omega_zeta(I_yy,C_M_a,C_M_da)
    % INPUT PARAMETERS
    % I_yy = 0.0075833;
    % C_M_a = -2e-5;
    % C_M_da = -1e-9;

    %C_M_da = input('Enter C_{M_{dot{alpha}}}: ');
    %C_M_a  = input('Enter C_{M_{alpha}}: ');
    %I_yy   = input('Enter I_{yy}: ');

    % Convert to standard 2nd-order form
    % ddot(alpha) + 2*zeta*omega_0*dot(alpha) + omega_0^2*alpha = 0
    a = -C_M_da / I_yy;
    b = -C_M_a  / I_yy;

    % Ensure system is physically meaningful
    if b < 0
        error('System is unstable (omega_0 imaginary). C_{M_alpha}/I_yy must be negative!');
    end

    omega_0 = sqrt(b);
    zeta = a / (2 * omega_0);
end

function plot_solution_at_angle(angle)
    ol_parameters = evalin('base', 'ol_parameters');
    control_surface_angles__rad = evalin('base', 'control_surface_angles__rad');
    zeta = interp1(control_surface_angles__rad , squeeze(ol_parameters(:,1,2)), angle);
    omega_0 = interp1(control_surface_angles__rad, squeeze(ol_parameters(:,1,1)), angle);

    omega_d = omega_0 * sqrt(1 - zeta^2);
    a_0 = 0.1;
    t = linspace(0, 10000000, 10000);
    alpha = a_0 * exp(-zeta * omega_0 * t) .* cos(omega_d * t);
    envelope = a_0 * exp(-zeta * omega_0 * t);
    figure;
    plot(t, alpha, 'b', 'LineWidth', 1.5);
    hold on;
    plot(t, envelope, 'r--', 'LineWidth', 1.5);
    plot(t, -envelope, 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('Angle [rad]');
    title(sprintf('Solution at angle %.2f rad', angle));
    legend('α(t)', 'Envelope +', 'Envelope -');
    grid on;

end
 
