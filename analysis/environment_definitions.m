
% Constants
altitude__m = 3e5;
gravitational_parameter__m3_per_s2 = 3.986e14;
radius__m = 6.378e6;

rotational_velocity_BI_B__rad_per_s = 0;
velocity_I_I__m_per_s = 7800 * [1;0;0];
wind_velocity_I_I__m_per_s = zeros(3,1);
%[T, R] = atmosnrlmsise00(altitude__m, 0, 0, 2024, 150, 0); -> Aerospace
%toolbox
density__kg_per_m3 = 3.8e-12;% R(6);
temperature__K = 934;%T(2);
particles_mass__kg = 16 * 1.6605390689252e-27;
temperature_ratio_method = 1;