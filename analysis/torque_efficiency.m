import vleo_aerodynamics_core.*
clear;

%% Setup environment, geometry, and LUT
geometry_type = 'shuttlecock'; % Options: 'plate' or 'shuttlecock'
lut_file = "aerodynamic_coefficients_panel_method_poly.mat";
[bodies, num_bodies, rotation_face_index, x_label, environment_definitions, lut_data] = setup(lut_file, geometry_type);
% Override rotation_face_index for this specific analysis
if strcmp(geometry_type, 'shuttlecock')
    rotation_face_index = 3;
end

%%
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles, 2);
aerodynamic_torque_B_B__Nm = aerodynamic_force_B__N;
aerodynamic_forces_B__N = nan(3,num_bodies,num_angles,2);
aerodynmaic_torques_B_B__Nm = aerodynamic_forces_B__N;
attitude_quaternion_BI = [1; 0; 0; 0];
for i = 1:num_angles
    for model = 1:2
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        [aerodynamic_force_B__N(:,i,model), ...
        aerodynamic_torque_B_B__Nm(:,i,model),...
        aerodynamic_forces_B__N(:,:,i,model),...
        aerodynmaic_torques_B_B__Nm(:,:,i,model)] = ...
        vleoAerodynamics(attitude_quaternion_BI, ...
                            environment_definitions.rotational_velocity_BI_B__rad_per_s, ...
                            environment_definitions.velocity_I_I__m_per_s, ...
                            environment_definitions.wind_velocity_I_I__m_per_s, ...
                            environment_definitions.density__kg_per_m3, ...
                            environment_definitions.temperature__K, ... 
                            environment_definitions.particles_mass__kg, ...
                            bodies, ...                                                       
                            bodies_rotation_angles__rad, ...
                            environment_definitions.temperature_ratio_method,...
                            model,...
                            lut_data);
    end
end

%% plot torque efficiency
figure;
hold on;
grid on;
plot(aerodynamic_torque_B_B__Nm(2,:,1), -aerodynamic_force_B__N(1,:,1), 'b', 'DisplayName', sprintf('Sentman \\alpha_E = %.4f',environment_definitions.energy_accommodation));
plot(aerodynamic_torque_B_B__Nm(2,:,2), -aerodynamic_force_B__N(1,:,2), 'r', 'DisplayName', 'IRS model');
xlabel('pitch Torque [Nm]');
ylabel('Drag [N]');
legend('Location', 'northwest');

%% Save figure as PNG and EPS
matlab2tikz(sprintf('torque_efficiency_%s.tex', geometry_type));

