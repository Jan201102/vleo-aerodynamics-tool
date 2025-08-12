%% LIFT TO DRAG RATIO
% calculate the lift to drag ratio for sentman and the new IRS model
% for a two sided flat plate and plot the results
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;
% Import constants from environment:definitions.m
run('environment_definitions.m');
energy_accommodation = 1;
temperature_ratio_method = 1;

%% load lut data
%lut_data = load_lut("aerodynamic_coefficients_panel_method.csv");
%lut_data = load("aerodynamic_coefficients_panel_method_sentman_spline.mat")
lut_data = load("aerodynamic_coefficients_panel_method_poly.mat");
%% load geometry based on parameter
bodies = parametrized_flat_plate(1, 1, [0,0,0],false,energy_accommodation,surface_temperature__K);
showBodies(bodies, [0], 0.75, 0.25);
num_bodies = 1;
rotation_face_index = 1;
x_label = "angle of attack [°]";

%% calculate aerodynamic forces

%loops

num_angles = 201;
control_surface_angles__rad = linspace(-pi/2, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        disp(current_angle)
        attitude_quarternion_BI = [[0;-1;0]*sin(current_angle/2); cos(current_angle/2)];
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        [aerodynamic_force_B__N(:,i,model), ~,~,~] = ...
            vleoAerodynamics(...
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
end

aerodynamic_force_sentman__N = aerodynamic_force_B__N(:,:,1);
aerodynamic_force_new__N = aerodynamic_force_B__N(:,:,2);

%% plot force envelopes for both models
figure;
hold on;
grid on;
plot3(aerodynamic_force_B__N(1,:,1), aerodynamic_force_B__N(2,:,1), aerodynamic_force_B__N(3,:,1), 'b', 'DisplayName', sprintf('Sentman Model \\alpha_E = %.4f',energy_accommodation));
plot3(aerodynamic_force_B__N(1,:,2), aerodynamic_force_B__N(2,:,2), aerodynamic_force_B__N(3,:,2), 'r', 'DisplayName', 'New IRS Model');
xlabel('X Force [N]');
ylabel('Y Force [N]');
zlabel('Z Force [N]');
title("Aerodynamic Force Envelopes for");
view([0 -1 0])
legend;

%% plot s force
figure;
hold on;
grid on;
plot(control_surface_angles__rad,aerodynamic_force_B__N(1,:,1), 'b', 'DisplayName',sprintf('Sentman Model \\alpha_E = %.4f',energy_accommodation))
plot(control_surface_angles__rad,aerodynamic_force_B__N(1,:,2), 'r', 'DisplayName', 'New IRS Model')
xlabel("angle")
ylabel("xforce")
legend()
%% Plot X Force vs Angle of Attack
figure;
hold on;
grid on;
control_surface_angles_plot = control_surface_angles__rad;
control_surface_angles_plot(control_surface_angles_plot > pi/4) = pi/2 - control_surface_angles_plot(control_surface_angles_plot > pi/4) ;
plot(control_surface_angles_plot,aerodynamic_force_B__N(1,:,1), 'b', 'DisplayName',sprintf('Sentman Model \\alpha_E = %.4f',energy_accommodation))
plot(control_surface_angles_plot,aerodynamic_force_B__N(1,:,2), 'r', 'DisplayName', 'New IRS Model')
xlabel("angle")
ylabel("xforce")
legend()

%% symmetry residuals for angle larger than 0, with symmetry at pi/4, 45°
x_sym = pi/4;
control_surface_angles_sym = control_surface_angles__rad(control_surface_angles__rad >= 0);
aerodynamic_force_filtered_B__N = aerodynamic_force_B__N(:,control_surface_angles__rad>=0,:);
plot_symmetry_residuals(control_surface_angles_sym, aerodynamic_force_filtered_B__N, x_sym);

%matlab2tikz('onesidedplate_symetry.tex')
%% symmetry towards 0
x_sym = 0;
plot_symmetry_residuals(control_surface_angles__rad,aerodynamic_force_B__N,x_sym)
%% Helper functions
function plot_symmetry_residuals(x, aerodynamic_force_B__N,x_sym)
    % PLOT_SYMMETRY_RESIDUALS Plot symmetry residuals for both Sentman and IRS models
    %
    % Inputs:
    %   x- angle array for plotting
    %   aerodynamic_force_B__N - force data (3 x num_angles x 2)
    %   x_sym - symmetry angle
    x_common = linspace(min(x), max(x), 100);
    
    % Get Y data for both models
    y_sentman = aerodynamic_force_B__N(1, :, 1);
    y_irs = aerodynamic_force_B__N(1, :, 2);

    % Interpolate both models to common x grid
    y_sentman_interp = interp1(x, y_sentman, x_common, 'linear', 'extrap');
    y_irs_interp = interp1(x, y_irs, x_common, 'linear', 'extrap');
    
    % Create mirrored versions
    x_mirrored = 2*x_sym - x_common;
    y_sentman_mirrored_interp = interp1(x, y_sentman, x_mirrored, 'linear', 'extrap');
    y_irs_mirrored_interp = interp1(x, y_irs, x_mirrored, 'linear', 'extrap');
    
    % Calculate residuals for both models
    residuals_sentman = y_sentman_interp - y_sentman_mirrored_interp;
    residuals_irs = y_irs_interp - y_irs_mirrored_interp;
    
    % Create the plot
    figure;
    subplot(2,1,1);
    plot(x_common, y_sentman_interp, 'b-');
    hold on;
    plot(x_common, y_sentman_mirrored_interp, 'b--');
    plot(x_common, y_irs_interp, 'r-');
    plot(x_common, y_irs_mirrored_interp, 'r--');
    xline(x_sym, 'k:', 'LineWidth', 1.5);
    legend('Sentman Original', 'Sentman Mirrored', 'IRS Original', 'IRS Mirrored', 'Symmetry Line');
    title('tangential forces');
    ylabel('tangential Force [N]');
    xlabel('Angle [°]');
    grid on;
    
    subplot(2,1,2);
    plot(x_common, residuals_sentman, 'b-', 'DisplayName','Sentman Model');
    hold on;
    plot(x_common, residuals_irs, 'r-', 'DisplayName', 'New IRS Model');
    xlabel('Angle [°]');
    ylabel('Residuals [N]');
    title('Symmetry Residuals (Original - Mirrored)');
    legend("Location","best");
    grid on;
end
