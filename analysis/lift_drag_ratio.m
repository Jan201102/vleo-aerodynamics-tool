%% LIFT TO DRAG RATIO
% calculate the lift to drag ratio for sentman and the new IRS model
% for a two sided flat plate and plot the results
import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;
run('environment_definitions.m');

%% Geometry selection parameter
geometry_type = 'plate'; % Options: 'plate' or 'shuttlecock'
temperature_ratio_method = 1;

%% load lut data
lut_data = load("aerodynamic_coefficients_panel_method_poly.mat");

%% load geometry based on parameter
[bodies, num_bodies, rotation_face_index, x_label] = load_geometry(geometry_type, energy_accommodation, surface_temperature__K);
if strcmp(geometry_type, 'shuttlecock')
    rotation_face_index = 3;
end

%% calculate aerodynamic forces
%loops
attitude_quarternion_BI = [1;0;0;0];
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
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

% calculate lift to drag ratio
lift_sentman__N = aerodynamic_force_sentman__N(3,:);
drag_sentman__N = -aerodynamic_force_sentman__N(1,:);
lift_new__N = aerodynamic_force_new__N(3,:);
drag_new__N = -aerodynamic_force_new__N(1,:);
lift_to_drag_ratio_sentman = lift_sentman__N./drag_sentman__N;
lift_to_drag_ratio_new = lift_new__N./drag_new__N;

%% plot results
figure;
hold on;
grid on;
plot(rad2deg(control_surface_angles__rad), lift_to_drag_ratio_sentman,"b", 'DisplayName', sprintf('Sentman Model \\alpha_E = %.4f',energy_accommodation));
plot(rad2deg(control_surface_angles__rad), lift_to_drag_ratio_new,"r", 'DisplayName', 'New IRS Model');
xlabel(x_label);
ylabel('lift to drag ratio');
legend;

%% Save figure as PNG and EPS
matlab2tikz(sprintf('lift_drag_ratio_%s.tex', geometry_type));

%% plot  drag an lift coefficients for both models in a 2x1 subplot grid
factor = 0.5*density__kg_per_m3*norm(velocity_I_I__m_per_s)^2;
figure;
subplot(2,1,1);
hold on;
grid on;
plot(rad2deg(control_surface_angles__rad), lift_sentman__N/factor, "b", 'DisplayName','Sentman Model');
plot(rad2deg(control_surface_angles__rad), lift_new__N/factor, "r", 'DisplayName', 'Schütte Model');
xlabel(x_label);
ylabel('C_L');
legend;
subplot(2,1,2);
hold on;
grid on;
plot(rad2deg(control_surface_angles__rad), drag_sentman__N/factor, "b", 'DisplayName',"Sentman Model");
plot(rad2deg(control_surface_angles__rad), drag_new__N/factor, "r", 'DisplayName', 'Schütte Model');
xlabel(x_label);
ylabel('C_D');
legend("location","northwest");
matlab2tikz('drag_lift_coeffcients.tikz')
