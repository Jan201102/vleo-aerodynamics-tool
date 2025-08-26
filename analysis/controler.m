%% CONTROLLER DESIGN AND ANALYSIS
% Design and analyze a pointing controller for a satellite with aerodynamic
% torques, comparing Sentman and IRS models across different densities
clear;
clc;

%% Load Aerodynamic Data
load('aerodynamic_stiffness_shuttlecock.mat');
load('aerodynamic_damping_shuttlecock.mat');

%% System Parameters
control_surface_angle = deg2rad(45);  % Control surface deflection angle [rad]
I_yy = 0.0375;                       % Moment of inertia about y-axis [kg*m2]

% Interpolate aerodynamic coefficients at control surface angle
c_m_sentman = interp1(control_surface_angles__rad, aero_stiffness(:,1), control_surface_angle);
c_m_irs = interp1(control_surface_angles__rad, aero_stiffness(:,2), control_surface_angle);
c_q_irs = interp1(control_surface_angles__rad, aero_damping(:,2), control_surface_angle);
c_q_sentman = interp1(control_surface_angles__rad, aero_damping(:,1), control_surface_angle);

%% Controller Design
T_epsilon = 130;                     % Settling time requirement [s]
zeta = 0.7;                         % Damping ratio
omega_0 = -log(0.02*sqrt(1-zeta^2))/(zeta*T_epsilon);  % Natural frequency [rad/s]
omega_0_deg = rad2deg(omega_0);     % Natural frequency [deg/s]

% Controller gains
k_p = I_yy * omega_0^2 + c_m_sentman;  % Proportional gain
k_d_min = 2 * I_yy * zeta * omega_0;   % Derivative gain

%% Transfer Functions and Closed-Loop Systems
G_sentman = tf([1], [I_yy, 0, -c_m_sentman]);
G_irs = tf([1], [I_yy, 0, -c_m_irs]);
K = tf([k_d_min, k_p], 1);

G_0_sentman = K * G_sentman;
G_0_irs = K * G_irs;
T_sentman = feedback(G_0_sentman, 1);
T_irs = feedback(G_0_irs, 1);

%% Bode Plot - Nominal Case
figure;
hold on;
bodeplot(G_0_sentman, 'b');
bodeplot(G_0_irs, 'r');
grid on;
legend('Sentman Model', 'Schütte Model');
title('Open-Loop Frequency Response');

% Display stability margins
[~, ~] = margin(G_0_sentman);
[~, ~] = margin(G_0_irs);
matlab2tikz('G0_bode_pointing.tikz', 'showInfo', false);

%% State-Space Representation
A_sentman = (1/I_yy) * [0, I_yy; -k_p+c_m_sentman, -k_d_min];
A_irs = (1/I_yy) * [0, I_yy; -k_p+c_m_irs, -k_d_min];
B = (1/I_yy) * [0, 0; k_p, k_d_min];
C = [1, 0];
D = 0;

% Initial conditions
initial_state_stable = [deg2rad(5), 0];     % 5° angular displacement
initial_state_rotating = [0, deg2rad(1)];   % 1°/s angular velocity

ss_sentman = ss(A_sentman, B, C, D);
ss_irs = ss(A_irs, B, C, D);

%% Initial Condition Response - Nominal Case
figure;

subplot(2,1,1);
[y, t] = initial(ss_sentman, initial_state_stable);
plot(t, rad2deg(y), 'b');
hold on;
grid on;
[y, t] = initial(ss_irs, initial_state_stable);
plot(t, rad2deg(y), 'r');
ylabel('Pitch Angle [°]');
xlabel('Time [s]');
xlim([0, max(t)]);
legend('Sentman Model', 'Schütte Model');
title('Stabilized - 5° Initial Displacement');

subplot(2,1,2);
[y, t] = initial(ss_sentman, initial_state_rotating);
plot(t, rad2deg(y), 'b');
hold on;
grid on;
[y, t] = initial(ss_irs, initial_state_rotating);
plot(t, rad2deg(y), 'r');
ylabel('Pitch Angle [°]');
xlabel('Time [s]');
xlim([0, max(t)]);
title('Rotating - 1°/s Initial Angular Velocity');
legend('Sentman Model', 'Schütte Model');
matlab2tikz('pointing_initial_condition_response.tikz');

%% Density Influence Study
particles_mass__kg = 16 * 1.6605390689252e-27;  % Atomic oxygen mass [kg]
n = 4.698e14;                                    % Number density [1/m3]
nominal_density__kg_per_m3 = particles_mass__kg * n;
density_max = 1.57e-10;                         % Maximum density [kg/m3]
density_min = 1.93e-13;                         % Minimum density [kg/m3]

% Scale aerodynamic coefficients with density
c_m_sentman_max = density_max / nominal_density__kg_per_m3 * c_m_sentman;
c_m_sentman_min = density_min / nominal_density__kg_per_m3 * c_m_sentman;
c_m_irs_max = density_max / nominal_density__kg_per_m3 * c_m_irs;
c_m_irs_min = density_min / nominal_density__kg_per_m3 * c_m_irs;

% Transfer functions for different densities
G_sentman_min = tf(1, [I_yy, 0, -c_m_sentman_min]);
G_irs_min = tf(1, [I_yy, 0, -c_m_irs_min]);
G_sentman_max = tf(1, [I_yy, 0, -c_m_sentman_max]);
G_irs_max = tf(1, [I_yy, 0, -c_m_irs_max]);

% Controller parameters for density extremes
T_epsilon = 130;
zeta = 0.707;
omega_0 = -log(0.02*sqrt(1-zeta^2)) / (zeta*T_epsilon);
k_p_max = 0;  % Proportional gain for maximum density
k_p_min = I_yy * omega_0^2 + c_m_sentman_min;
k_d_min = 2 * I_yy * zeta * omega_0;
k_d_max = 2 * I_yy * zeta * sqrt(-c_m_sentman_max/I_yy);

% Controllers and open-loop systems for density extremes
K_max = tf([k_d_max, k_p_max], 1);
K_min = tf([k_d_min, k_p_min], 1);
G_0_sentman_max = K_max * G_sentman_max;
G_0_sentman_min = K_min * G_sentman_min;
G_0_irs_max = K_max * G_irs_max;
G_0_irs_min = K_min * G_irs_min;
T_sentman_max = feedback(G_0_sentman_max, 1);
T_irs_max = feedback(G_0_irs_max, 1);
T_sentman_min = feedback(G_0_sentman_min, 1);
T_irs_min = feedback(G_0_irs_min, 1);

% Open-loop parameters
omega_0_sentman_max = sqrt(-c_m_sentman_max/I_yy);
omega_0_sentman_min = sqrt(-c_m_sentman_min/I_yy);
omega_0_irs_max = sqrt(-c_m_irs_max/I_yy);
omega_0_irs_min = sqrt(-c_m_irs_min/I_yy);

%% State-Space Representation for Density Study
A_sentman_max = (1/I_yy) * [0, I_yy; -k_p_max+c_m_sentman_max, -k_d_max];
A_irs_max = (1/I_yy) * [0, I_yy; -k_p_max+c_m_irs_max, -k_d_max];
A_sentman_min = (1/I_yy) * [0, I_yy; -k_p_min+c_m_sentman_min, -k_d_min];
A_irs_min = (1/I_yy) * [0, I_yy; -k_p_min+c_m_irs_min, -k_d_min];
B_max = (1/I_yy) * [0, 0; k_p_max, k_d_max];
B_min = (1/I_yy) * [0, 0; k_p_min, k_d_min];
C = [1, 0];
D = 0;

% Initial conditions
initial_state_stable = [deg2rad(5), 0];
initial_state_rotating = [0, deg2rad(1)];

ss_sentman_max = ss(A_sentman_max, B_max, C, D);
ss_irs_max = ss(A_irs_max, B_max, C, D);
ss_sentman_min = ss(A_sentman_min, B_min, C, D);
ss_irs_min = ss(A_irs_min, B_min, C, D);
%% Bode Plot - Density Study
figure;
hold on;
bodeplot(G_0_sentman_min, 'b');
bodeplot(G_0_sentman_max, 'b--');
bodeplot(G_0_irs_min, 'r');
bodeplot(G_0_irs_max, 'r--');
grid on;

legend(sprintf('Sentman \\rho = %.2e', density_min), ...
       sprintf('Sentman \\rho = %.2e', density_max), ...
       sprintf('IRS \\rho = %.2e', density_min), ...
       sprintf('IRS \\rho = %.2e', density_max));
title('Open-Loop Frequency Response - Density Study');
matlab2tikz('G0_bode_density_study.tex', 'showInfo', false);

% Display stability margins for all cases
[~, ~] = margin(G_0_sentman_max);
[~, ~] = margin(G_0_sentman_min);
[~, ~] = margin(G_0_irs_max);
[~, ~] = margin(G_0_irs_min);

%% Initial Condition Response - Density Study
figure;

subplot(2,1,1);
hold on;
grid on;
[y,t] = initial(ss_sentman_min,initial_state_stable);
plot(t,rad2deg(y),'b')
[y,t] = initial(ss_sentman_max,initial_state_stable);
plot(t,rad2deg(y),'b--')
[y,t] = initial(ss_irs_min,initial_state_stable);
plot(t,rad2deg(y),'r')
xlim([0 max(t)])
[y,t] = initial(ss_irs_max,initial_state_stable);
plot(t,rad2deg(y),'r--')
ylabel('Pitch Angle [°]')
xlabel('Time [s]')

legend(sprintf('sentman \\rho = %.2e',density_min),...
    sprintf('sentman \\rho = %.2e',density_max),...
    sprintf('irs \\rho = %.2e',density_min),...
    sprintf('irs \\rho = %.2e',density_max));
title('Stabilized')

subplot(2,1,2);
hold on;
grid on;
[y,t] = initial(ss_sentman_min,initial_state_rotating);
plot(t,rad2deg(y),'b')
[y,t] = initial(ss_sentman_max,initial_state_rotating);
plot(t,rad2deg(y),'b--')
[y,t] = initial(ss_irs_min,initial_state_rotating);
plot(t,rad2deg(y),'r')
xlim([0 max(t)])
[y,t] = initial(ss_irs_max,initial_state_rotating);
plot(t,rad2deg(y),'r--')
ylabel('Pitch Angle [°]')
xlabel('Time [s]')

title('Rotating')
legend(sprintf('sentman \\rho = %.2e',density_min),...
    sprintf('sentman \\rho = %.2e',density_max),...
    sprintf('irs \\rho = %.2e',density_min),...
    sprintf('irs \\rho = %.2e',density_max));

matlab2tikz('pointing_initial_condition_response_density.tikz');
