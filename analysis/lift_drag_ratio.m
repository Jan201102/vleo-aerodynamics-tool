%% LIFT TO DRAG RATIO
% calculate the lift to drag ratio for sentman and the new IRS model
% for a two sided flat plate and plot the results
import vleo_aerodynamics_core.*
clear;

%% Setup environment, geometry, and LUT
geometry_type = 'plate'; % Options: 'plate' or 'shuttlecock'
lut_file = "aerodynamic_coefficients_panel_method_poly.mat";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);
if strcmp(geometry_type, 'shuttlecock')
    rotation_face_index = 3;
end

%% calculate aerodynamic forces
%loops
attitude_quaternion_BI = [1;0;0;0];
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
plot(rad2deg(control_surface_angles__rad), lift_to_drag_ratio_sentman,"b", 'DisplayName', sprintf('Sentman Model \\alpha_E = %.4f', environment_definitions.energy_accommodation));
plot(rad2deg(control_surface_angles__rad), lift_to_drag_ratio_new,"r", 'DisplayName', 'New IRS Model');
xlabel(x_label);
ylabel('lift to drag ratio');
legend;

%% Save figure as PNG and EPS
matlab2tikz(sprintf('lift_drag_ratio_%s.tex', geometry_type));

%% plot  drag an lift coefficients for both models in a 2x1 subplot grid
factor = 0.5*environment_definitions.density__kg_per_m3*norm(environment_definitions.velocity_I_I__m_per_s)^2;
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
