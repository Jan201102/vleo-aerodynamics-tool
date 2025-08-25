clear
clc

load("aerodynamic_stiffness_shuttlecock.mat")
load("aerodynamic_dampening_shuttlecock.mat")

control_surface_angle = deg2rad(45);
c_m_sentman = interp1(control_surface_angles__rad,aero_stiffness(:,1),control_surface_angle)
c_m_irs = interp1(control_surface_angles__rad,aero_stiffness(:,2),control_surface_angle)
c_q_irs = interp1(control_surface_angles__rad,aero_dampening(:,2),control_surface_angle)
c_q_sentman = interp1(control_surface_angles__rad,aero_dampening(:,1),control_surface_angle)
I_yy = 0.0375;

%% controler
T_epsilon = 130;
zeta = 0.7;
omega_0 =-log(0.02*sqrt(1-zeta^2))/(zeta*T_epsilon)
omega_0_deg = rad2deg(omega_0)
k_p = I_yy * omega_0^2 + c_m_sentman
k_d_min = 2* I_yy * zeta*omega_0

G_sentman = tf([1],[I_yy,0,-c_m_sentman]);
G_irs = tf([1],[I_yy,0,-c_m_irs]);
K = tf([k_d_min,k_p],1)
G_0_sentman = K*G_sentman;
G_0_irs = K*G_irs;
T_sentman = feedback(G_0_sentman,1);%tf([k_d,k_p],[I_yy,k_d, (k_p-c_m_sentman)]);
T_irs = feedback(G_0_irs,1);%tf([k_d,k_p],[I_yy,k_d, (k_p-c_m_irs)]);

%% bode plot
figure
hold on
bodeplot(G_0_sentman,"b");
bodeplot(G_0_irs,"r");
grid on;
legend("Sentman Model","Schütte Model")
[g,p] = margin(G_0_sentman)
[g,p] = margin(G_0_irs)
matlab2tikz('G0_bode_pointing.tikz','showInfo', false);

%% response to inital conditions
A_sentman = 1/I_yy.* [0 I_yy; -k_p+c_m_sentman -k_d_min];
A_irs = 1/I_yy.* [0 I_yy; -k_p+c_m_irs -k_d_min];
B = 1/I_yy.*[0 0; k_p k_d_min];
C = [1 0];
D = 0;
initial_state_stable = [ deg2rad(5) 0];
initial_state_rotating = [0 deg2rad(1)];
ss_sentman = ss(A_sentman,B,C,D);
ss_irs = ss(A_irs,B,C,D);

figure
subplot(2,1,1)
[y,t] = initial(ss_sentman,initial_state_stable);
plot(t,rad2deg(y),"b")
hold on;
grid on;
[y,t] = initial(ss_irs,initial_state_stable);
plot(t,rad2deg(y),"r")
ylabel("Pitch Angle [°]")
xlabel("Time [s]")
xlim([0 max(t)])
legend("Sentman Model","Schütte Model")
title("Stabilized")
subplot(2,1,2)
[y,t] = initial(ss_sentman,initial_state_rotating);
plot(t,rad2deg(y),"b")
hold on;
grid on;
[y,t] = initial(ss_irs,initial_state_rotating);
plot(t,rad2deg(y),"r")
ylabel("Pitch Angle [°]")
xlabel("Time [s]")
xlim([0 max(t)])
title("Rotating")
legend("Sentman Model","Schütte Model")
matlab2tikz('pointing_initial_condition_response.tikz');

%% Density influence
particles_mass__kg = 16 * 1.6605390689252e-27;
n = 4.698e14;
nominal_density__kg_per_m3 = particles_mass__kg * n;
density_max = 1.57e-10; %particles_mass__kg*n_max
density_min = 1.93e-13; %particles_mass__kg*n_min

c_m_sentman_max = density_max/nominal_density__kg_per_m3 * c_m_sentman;
c_m_sentman_min = density_min/nominal_density__kg_per_m3 * c_m_sentman;
c_m_irs_max = density_max/nominal_density__kg_per_m3 * c_m_irs;
c_m_irs_min = density_min/nominal_density__kg_per_m3 *c_m_irs;

G_sentman_min = tf(1,[I_yy,0,-c_m_sentman_min]);
G_irs_min = tf(1,[I_yy,0,-c_m_irs_min]);
G_sentman_max = tf(1,[I_yy,0,-c_m_sentman_max]);
G_irs_max = tf(1,[I_yy,0,-c_m_irs_max]);

T_epsilon = 130;
zeta = 0.707;
omega_0 =-log(0.02*sqrt(1-zeta^2))/(zeta*T_epsilon);
omega_0_deg = rad2deg(omega_0);
k_p_max =0;%I_yy * omega_0^2 + c_m_sentman_max;
k_p_min = I_yy *omega_0^2 +c_m_sentman_min;
k_d_min = 2* I_yy * zeta*omega_0;
k_d_max = 2*I_yy*zeta*sqrt(-c_m_sentman_max/I_yy);;

K_max = tf([k_d_max,k_p_max],1);
K_min = tf([k_d_min,k_p_min],1);
G_0_sentman_max = K_max*G_sentman_max;
G_0_sentman_min = K_min*G_sentman_min;
G_0_irs_max = K_max*G_irs_max;
G_0_irs_min = K_min*G_irs_min;
T_sentman_max = feedback(G_0_sentman_max,1);
T_irs_max = feedback(G_0_irs_max,1);
T_sentman_min = feedback(G_0_sentman_min,1);
T_irs_min = feedback(G_0_irs_min,1);

%Ol parameter
omega_0_sentman_max = sqrt(-c_m_sentman_max/I_yy);
omega_0_sentman_min = sqrt(-c_m_sentman_min/I_yy);
omega_0_irs_max = sqrt(-c_m_irs_max/I_yy);
omega_0_irs_min = sqrt(-c_m_irs_min/I_yy);

%state-space representation of closed loop for different models and densities
A_sentman_max = 1/I_yy.* [0 I_yy; -k_p_max+c_m_sentman_max -k_d_max];
A_irs_max = 1/I_yy.* [0 I_yy; -k_p_max+c_m_irs_max -k_d_max];
A_sentman_min = 1/I_yy.* [0 I_yy; -k_p_min+c_m_sentman_min -k_d_min];
A_irs_min = 1/I_yy.* [0 I_yy; -k_p_min+c_m_irs_min -k_d_min];
B_max = 1/I_yy.*[0 0; k_p_max k_d_max];
B_min = 1/I_yy.*[0 0; k_p_min k_d_min];
C = [1 0];
D = 0;
initial_state_stable = [ deg2rad(5) 0];
initial_state_rotating = [0 deg2rad(1)];
ss_sentman_max = ss(A_sentman_max,B_max,C,D);
ss_irs_max = ss(A_irs_max,B_max,C,D);
ss_sentman_min = ss(A_sentman_min,B_min,C,D);
ss_irs_min = ss(A_irs_min,B_min,C,D);
%% bodeplot
figure
hold on
bodeplot(G_0_sentman_min,"b");
bodeplot(G_0_sentman_max,"b--");
bodeplot(G_0_irs_min,"r");
bodeplot(G_0_irs_max,"r--");
grid on;

legend(sprintf("sentman \\rho = %.2e",density_min),...
    sprintf("sentman \\rho = %.2e",density_max),...
    sprintf("irs \\rho = %.2e",density_min),...
    sprintf("irs \\rho = %.2e",density_max))
matlab2tikz('G0_bode_density_study.tex','showInfo', false);


[~,~] = margin(G_0_sentman_max)
[~,~] = margin(G_0_sentman_min)
[~,~] = margin(G_0_irs_max)
[~,~] = margin(G_0_irs_min)

%% response to initial conditions
figure
subplot(2,1,1)
hold on;
grid on;
[y,t] = initial(ss_sentman_min,initial_state_stable);
plot(t,rad2deg(y),"b")
[y,t] = initial(ss_sentman_max,initial_state_stable);
plot(t,rad2deg(y),"b--")
[y,t] = initial(ss_irs_min,initial_state_stable);
plot(t,rad2deg(y),"r")
xlim([0 max(t)])
[y,t] = initial(ss_irs_max,initial_state_stable);
plot(t,rad2deg(y),"r--")
ylabel("Pitch Angle [°]")
xlabel("Time [s]")

legend(sprintf("sentman \\rho = %.2e",density_min),...
    sprintf("sentman \\rho = %.2e",density_max),...
    sprintf("irs \\rho = %.2e",density_min),...
    sprintf("irs \\rho = %.2e",density_max));
title("Stabilized")

subplot(2,1,2)
hold on;
grid on;
[y,t] = initial(ss_sentman_min,initial_state_rotating);
plot(t,rad2deg(y),"b")
[y,t] = initial(ss_sentman_max,initial_state_rotating);
plot(t,rad2deg(y),"b--")
[y,t] = initial(ss_irs_min,initial_state_rotating);
plot(t,rad2deg(y),"r")
xlim([0 max(t)])
[y,t] = initial(ss_irs_max,initial_state_rotating);
plot(t,rad2deg(y),"r--")
ylabel("Pitch Angle [°]")
xlabel("Time [s]")

title("Rotating")
legend(sprintf("sentman \\rho = %.2e",density_min),...
    sprintf("sentman \\rho = %.2e",density_max),...
    sprintf("irs \\rho = %.2e",density_min),...
    sprintf("irs \\rho = %.2e",density_max));

matlab2tikz('pointing_initial_condition_response_density.tikz');
