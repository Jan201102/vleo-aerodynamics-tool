%% CUBE MOMENTS ANALYSIS
% Calculate aerodynamic moments for cube geometry across different attitudes
% Compare Sentman and IRS models for attitude-dependent behavior
import vleo_aerodynamics_core.*
clear;

%% Setup Environment, Geometry, and LUT  
geometry_type = 'plate';
lut_file = 'aerodynamic_coefficients_panel_method_poly.mat';
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);

%% Calculate Aerodynamic Forces Across Attitudes
num_angles = 201;
control_surface_angles__rad = linspace(-pi/2, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles, 2);

for model = 1:2
    fprintf('Processing model %d...\n', model);
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        
        % Create attitude quaternion for rotation about y-axis
        attitude_quaternion_BI = [cos(current_angle/2); 0; -sin(current_angle/2); 0];
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        
        [aerodynamic_force_B__N(:,i,model), ~, ~, ~] = ...
            vleoAerodynamics(...
                attitude_quaternion_BI, ...
                environment_definitions.rotational_velocity_BI_B__rad_per_s, ...
                environment_definitions.velocity_I_I__m_per_s, ...
                environment_definitions.wind_velocity_I_I__m_per_s, ...
                environment_definitions.density__kg_per_m3, ...
                environment_definitions.temperature__K,...
                environment_definitions.particles_mass__kg,...
                bodies,...
                bodies_rotation_angles__rad,...
                environment_definitions.temperature_ratio_method,...
                model,...
                lut_data);
    end
end

aerodynamic_force_sentman__N = aerodynamic_force_B__N(:,:,1);
aerodynamic_force_new__N = aerodynamic_force_B__N(:,:,2);

%% symmetry residuals for angle larger than 0, with symmetry at pi/4, 45°
x_sym = pi/4;
control_surface_angles_sym = control_surface_angles__rad(control_surface_angles__rad >= 0);
aerodynamic_force_filtered_B__N = aerodynamic_force_B__N(:,control_surface_angles__rad>=0,:);
plot_symmetry_residuals(control_surface_angles_sym, aerodynamic_force_filtered_B__N, x_sym);

%matlab2tikz('onesidedplate_symetry.tex')

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
    legend('Location','best');
    grid on;
end
