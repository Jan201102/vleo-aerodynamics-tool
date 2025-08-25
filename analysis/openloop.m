%% OPEN LOOP ANALYSIS
% calculate the aerodynamic stiffness and damping for open loop analysis
import vleo_aerodynamics_core.*
clear;

%% Setup environment, geometry, and LUT
geometry_type = 'shuttlecock';
lut_file = "aerodynamic_coefficients_panel_method_poly.mat";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);
I_yy = 0.0375; % Moment of inertia around y-axis

%% Derivation Configuration
derivation_method = 'central';  % Options: 'central', 'five_point_stencil', 'seven_point_stencil', 'forward', 'backward'
delta_q_rad_per_s = 1e2;                  % Step size for numerical differentiation
delta_alpha = 1e-4;                  % Step size for numerical differentiation
axis_direction = [0; 1; 0];               % Rotation axis (y-axis for pitch)
torque_component = 2;                      % Component index (2 = y-component for pitch)

%% calculate aerodynamic stiffness and damping
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
ol_parameters = zeros(num_angles,2,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        
        C_M_da = calculate_aerodynamic_damping(derivation_method, ...
                                                delta_q_rad_per_s, ...
                                                axis_direction,...
                                                torque_component, ...
                                                bodies, ...
                                                bodies_rotation_angles__rad, ...
                                                model, ...
                                                lut_data, ...
                                                environment_definitions);
        C_M_a = calculate_aerodynamic_stiffness(derivation_method, ...
                                                0,...
                                                delta_alpha, ...
                                                axis_direction,...
                                                torque_component, ...
                                                bodies, ...
                                                bodies_rotation_angles__rad, ...
                                                model, ...
                                                lut_data, ...
                                                environment_definitions);

        omega_0 =  sqrt(-C_M_a/I_yy);
        ol_parameters(i,model,1) =omega_0;
        ol_parameters(i,model,2) = -C_M_da/(2*I_yy*omega_0);
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

%% compare omega_0 and omega_d by plotting relative difference over control surface angle
omega_d = ol_parameters(:,:,1).* sqrt(1 - ol_parameters(:,:,2).^2);
omega_0 = ol_parameters(:,:,1);
relative_difference = (omega_0 - omega_d) ./ omega_0;
figure;
plot(control_surface_angles__rad * 180 / pi, relative_difference(:,1), 'b-');
hold on;
plot(control_surface_angles__rad * 180 / pi, relative_difference(:,2), 'r-');
grid on;
xlabel(x_label);
ylabel('Relative Difference [-]');
title('Relative Difference between \omega_0 and \omega_d für verschiedene Modelle');
legend("Sentman","IRS","Location","southeast");


 
