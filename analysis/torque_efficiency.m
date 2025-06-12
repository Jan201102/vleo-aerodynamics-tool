import vleo_aerodynamics_core.*
addpath("analysis/functions");
clear;

% Import constants from environment:definitions.m
run('environment_definitions.m');

%% Geometry selection parameter
% Set to 'plate' for flat plate geometry or 'shuttlecock' for shuttlecock model
geometry_type = 'plate'; % Options: 'plate' or 'shuttlecock'
energy_accommodation = 0.0;
%% load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'aerodynamic_coefficients_panel_method.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%% load geometry based on parameter
if strcmp(geometry_type, 'plate')
    bodies = parametrized_flat_plate(0.13333333, 0.098, 0.07, 0.0,false,energy_accommodation);
    showBodies(bodies, [0], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
elseif strcmp(geometry_type, 'shuttlecock')
    bodies = load_from_gmsh(energy_accommodation);
    showBodies(bodies, [0,0/4,pi/4,pi/4,pi/4], 0.75, 0.25);
    num_bodies = 5;
    rotation_face_index = 4;
elseif strcmp(geometry_type, 'shuttlecock_wing')
    bodies = load_shuttlecock_wing();
    showBodies(bodies, [0/4], 0.75, 0.25);
    num_bodies = 1;
    rotation_face_index = 1;
else
    error("Invalid geometry_type. Use 'plate' or 'shuttlecock'.");
end

%%
num_angles = 101;
control_surface_angles__rad = linspace(0, pi/2, num_angles);
aerodynamic_force_B__N = nan(3, num_angles, 2);
aerodynamic_torque_B_B__Nm = aerodynamic_force_B__N;
attitude_quaternion_BI = [1; 0; 0; 0];
for i = 1:num_angles
    for model = 1:2
        current_angle = control_surface_angles__rad(i);
        bodies_rotation_angles__rad = zeros(1, num_bodies);
        bodies_rotation_angles__rad(rotation_face_index) = current_angle;
        [aerodynamic_force_B__N(:,i,model), ...
            aerodynamic_torque_B_B__Nm(:,i,model)] = ...
            vleoAerodynamics(attitude_quaternion_BI, ...
                             rotational_velocity_BI_B__rad_per_s, ...
                             velocity_I_I__m_per_s, ...
                             wind_velocity_I_I__m_per_s, ...
                             density__kg_per_m3, ...
                             temperature__K, ... 
                             particles_mass__kg, ...
                             bodies, ...                                                       
                             bodies_rotation_angles__rad, ...
                             temperature_ratio_method,...
                             model,...
                             2, ...
                             lut);
    end
end
%% plot force envelopes for each bodie in a subplot with tiled layout
if ndims(aerodynamic_torque_B_B__Nm) == 4
    figure;
    tl = tiledlayout('flow');
    title(tl, sprintf('Aerodynamic Forces for %s geometry', geometry_type));
    for b = 1:num_bodies
        ax = nexttile;
        grid on;
        hold on;
        title(ax, sprintf('Body %d', b));
        xlabel("x");
        ylabel("y");
        zlabel("z [N]");
        plot3(ax, squeeze(aerodynamic_force_B__N(1,b,:,1)), squeeze(aerodynamic_force_B__N(2,b,:,1)), squeeze(aerodynamic_force_B__N(3,b,:,1)), "b-o", 'DisplayName', 'Sentman');
        plot3(ax, squeeze(aerodynamic_force_B__N(1,b,:,2)), squeeze(aerodynamic_force_B__N(2,b,:,2)), squeeze(aerodynamic_force_B__N(3,b,:,2)), "r-o", 'DisplayName', 'IRS');
        legend(ax,"location","northwest");
        view(ax, [0 -1 0]);
    end
end

%% plot torque efficiency
figure;
title(sprintf('Aerodynamic Torque efficiency for %s geometry', geometry_type));
hold on;
grid on;
plot(aerodynamic_torque_B_B__Nm(2,:,1), -aerodynamic_force_B__N(1,:,1), 'b', 'DisplayName', sprintf('Sentman \\alpha_E = %.2f',energy_accommodation));
plot(aerodynamic_torque_B_B__Nm(2,:,2), -aerodynamic_force_B__N(1,:,2), 'r', 'DisplayName', 'IRS model');
xlabel('pitch Torque [Nm]');
ylabel('Drag [N]');
legend('Location', 'northwest');

%% Save figure as PNG and EPS
saveas(gcf, sprintf('torque_efficiency_%s.png', geometry_type));
saveas(gcf, sprintf('torque_efficiency_%s.eps',geometry_type), 'epsc');
%% plot force envelope

figure;
tl = tiledlayout('flow');
title(tl, sprintf('Aerodynamic Forces and Torques for %s gemoetry', geometry_type));
ax1 = nexttile;
grid on;
hold on;
title(ax1, 'Aerodynamic Forces');
xlabel("x");
ylabel("y");
zlabel("z [N]");
plot3(ax1,aerodynamic_force_B__N(1,:,1), aerodynamic_force_B__N(2,:,1), aerodynamic_force_B__N(3,:,1),"b",'DisplayName', 'Sentman');
plot3(ax1,aerodynamic_force_B__N(1,:,2), aerodynamic_force_B__N(2,:,2), aerodynamic_force_B__N(3,:,2),"r",'DisplayName', 'IRS');
legend(ax1);
view(ax1,[0 -1 0]);

ax2 = nexttile;
grid on;
hold on;
title(ax2, 'Aerodynamic Torques');
xlabel("x");
ylabel("y");
zlabel("z [Nm]");
plot3(ax2,aerodynamic_torque_B_B__Nm(1,:,1), aerodynamic_torque_B_B__Nm(2,:,1), aerodynamic_torque_B_B__Nm(3,:,1),"b",'DisplayName', 'Sentman');
plot3(ax2,aerodynamic_torque_B_B__Nm(1,:,2), aerodynamic_torque_B_B__Nm(2,:,2), aerodynamic_torque_B_B__Nm(3,:,2),"r",'DisplayName', 'IRS');
legend(ax2);
view(ax2,[0 -1 0]);

