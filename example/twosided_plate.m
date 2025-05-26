import vleo_aerodynamics_core.*
clear;


%load model data
[test_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
display(test_folder)
lut = fullfile(test_folder, 'cl_cd_cVAE_A01_flat_and_bird.csv');
if ~isfile(lut)
    error("Look-up table file not found. Please check the path: %s", lut);
end

%create a unit plate of two triangles laying in the xy plane
%with the normal pointing in the negative z direction
% and the centroid at the origin
plate = cell(1,1);
plate{1}.vertices_B = zeros(3,3,4);
plate{1}.vertices_B(:,:,1) = [-0.5 0.5 -0.5;
                            -0.5 -0.5 0.5;
                            -1e-10 -1e-10 -1e-10];
plate{1}.vertices_B(:,:,2) = [0.5 0.5 -0.5;
                            0.5 -0.5 0.5;
                            -1e-10 -1e-10 -1e-10];


 plate{1}.vertices_B(:,:,3) = [-0.5 0.5 -0.5;
                            -0.5 -0.5 0.5;
                            1e-10 1e-10 1e-10];
plate{1}.vertices_B(:,:,4) = [0.5 0.5 -0.5;
                            0.5 -0.5 0.5;
                            1e-10 1e-10 1e-10]; 



plate{1}.centroids_B = squeeze(mean(plate{1}.vertices_B, 2));
plate{1}.normals_B = [0 0 0 0; 0 0 0 0; -1 -1 1 1];
plate{1}.areas = [0.5 0.5 0.5 0.5];
%rotation hinge point is the centroid of the plate
plate{1}.rotation_hinge_point_B = [0;0;0];
plate{1}.rotation_direction_B = [0;-1;0];
plate{1}.temperatures__K = [300;300;300;300];
plate{1}.energy_accommodation_coefficients = [0.9;0.9;0.9;0.9];

%visualize body
showBodies(plate, [pi/4], 0.75, 0.25);

%set parameters for orbit,atmosphere and methods
altitude__m = 3e5;
gravitational_parameter__m3_per_s2 = 3.986e14;
radius__m = 6.378e6;

rotational_velocity_BI_B__rad_per_s = 0;
velocity_I_I__m_per_s = sqrt(gravitational_parameter__m3_per_s2 ...
                             / (radius__m + altitude__m)) * [1;0;0];
wind_velocity_I_I__m_per_s = zeros(3,1);
%[T, R] = atmosnrlmsise00(altitude__m, 0, 0, 2024, 150, 0); -> Aerospace toolbox
density__kg_per_m3 = 3.8e-12;% R(6);
temperature__K = 950;%T(2);
particles_mass__kg = 16 * 1.6605390689252e-27;
temperature_ratio_method = 1;

%loops
attitude_quarternion_BI = [1;0;0;0];
num_angles = 101;
control_surface_angles__rad = linspace(0, pi, num_angles);
aerodynamic_force_B__N = nan(3, num_angles,2);
for model = 1:2
    for i = 1:num_angles
        current_angle = control_surface_angles__rad(i);
        [aerodynamic_force_B__N(:,i,model), ~] = ...
            vleoAerodynamics(...
                attitude_quarternion_BI,...
                rotational_velocity_BI_B__rad_per_s,...
                velocity_I_I__m_per_s,...
                wind_velocity_I_I__m_per_s,...
                density__kg_per_m3,...
                temperature__K,...
                particles_mass__kg,...
                plate,...
                current_angle,...
                temperature_ratio_method,...
                model,...
                lut);
    end
end
%plot forces
figure;
tl = tiledlayout('flow');
title(tl, 'Aerodynamic Forces and Torques for a two-sided plate');
ax1 = nexttile;
grid on;
hold on;
title(ax1, 'Aerodynamic Forces');
xlabel("x");
ylabel("y");
zlabel("z [N]");
plot3(ax1,aerodynamic_force_B__N(1,:,1), aerodynamic_force_B__N(2,:,1), aerodynamic_force_B__N(3,:,1),'DisplayName', 'Sentmann');
plot3(ax1,aerodynamic_force_B__N(1,:,2), aerodynamic_force_B__N(2,:,2), aerodynamic_force_B__N(3,:,2),'DisplayName', 'IRS');
legend(ax1);
view(ax1,[0 -1 0])
