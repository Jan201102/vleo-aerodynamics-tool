% Constants

%DSMC parameter:
% T_i=934K
% T_w = 300K
% u = 7800 m/s
% n = 4.698e14 1/m^3
% partikel_mass = 16U (Atomarer auserstoff) -> rho = 1.2482e-11 kg/m^3
rotational_velocity_BI_B__rad_per_s = 0;
velocity_I_I__m_per_s = 7800 * [1;0;0];
wind_velocity_I_I__m_per_s = zeros(3,1);
n = 4.698e14;
temperature__K = 934;
surface_temperature__K = 300;
particles_mass__kg = 16 * 1.6605390689252e-27;
density__kg_per_m3 = particles_mass__kg * n;
temperature_ratio_method = 1;
energy_accommodation = (7.5e-17 * n * temperature__K)/(1+7.5e-17 * n * temperature__K);