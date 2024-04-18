% Author: Minh Phat Do
% Description: Consensus-based CACC, constant graph

%% Reset workspace
clear;
clc;
close all;

%% System definition

% Define the plattoon configuration
G = digraph([1 2 3], [2 3 4]);
% plot(G)
Adj = adjacency(G)*eye(4);
S = [1;1;1;1];
L0 = Adj*S;
D0 = diag(L0);
L = D0-Adj;
P = diag([0 0 0 1]);
L_hat = L + P
lamda_L_hat = eig(L_hat);
lamda_2 = lamda_L_hat(2);

%% System definition
lv = 4.46;
h = 1;
rv = 0;
% q0_init = 0;
% q1_init = 0;
% q2_init = 0;
% q3_init = 0;
% q4_init = 0;
% q5_init = 0;
% q6_init = 0;
% q7_init = 0;
% q8_init = 0;
% q9_init = 0;
% q10_init = 0;
range = 10;
q0_init = rand() * range;
q1_init = rand() * range;
q2_init = rand() * range;
q3_init = rand() * range;
q4_init = rand() * range;
q5_init = rand() * range;
q6_init = rand() * range;
q7_init = rand() * range;
q8_init = rand() * range;
q9_init = rand() * range;
q10_init = rand() * range;

tau = 0.1;
A = [0 1 0; 0 0 1; 0 0 -1/tau];
B = [0 0 1/tau]';
C = diag([1 1 1]);
D = [0 0 0]';

error_ss = ss(A, B, C, D);
Q = diag([0.09 0.5 1]);
R = 1;
K_lqr = lqr(error_ss, Q, R);

%% LMIs setup

% Define SDP variables
X = sdpvar(3,3);
Y = sdpvar(1,3,'full');

% Define contraints and solve LMIs
kappa = 1;
F = [X>=0];
F = [F, X*A' + A*X + lamda_2*B*Y + conj(lamda_2)*Y'*B' + 2*kappa*X<= 0];
options = sdpsettings();
options.verbose = 0;
optimize(F,0,options);
K_lmi = value(Y)*inv(value(X));

% K_cacc = [0.2 1.2 0] %from the author

% K_cacc = [K_lqr(1) K_lqr(2) 0];
K_cacc = -K_lmi

%% Results

% Test on error model
% 3 vehicles
% sim("fixed_graph_sim.slx")
% 10 vehicles
sim("fixed_graph_sim_compact.slx")


sim_time = tout;

% Spacing error
e_1 = error.Data(:,1);
e_2 = error.Data(:,2);
e_3 = error.Data(:,3);
e_4 = error.Data(:,4);
e_5 = error.Data(:,5);
e_6 = error.Data(:,6);
e_7 = error.Data(:,7);
e_8 = error.Data(:,8);
e_9 = error.Data(:,9);
e_10 = error.Data(:,10);

lw = 0.5;

figure(1)
subplot(212)
plot(sim_time, e_1, ...
    sim_time, e_2, ...
    sim_time, e_3, ...
    sim_time, e_4, ...
    sim_time, e_5, ...
    sim_time, e_6, ...
    sim_time, e_7, ...
    sim_time, e_8, ...
    sim_time, e_9, ...
    sim_time, e_10, ...
    LineWidth=lw); xlim([0 20]); ylim([-15 15])
% legend("e1", "e2", "e3");
xlabel('t(s)')
ylabel('$e_i$(m)', 'Interpreter','latex')
grid on

% Position
% p_0 = monitor.Data(:,4);
% p_1 = monitor.Data(:,5);
% p_2 = monitor.Data(:,6);
% p_3 = monitor.Data(:,7);
% subplot(222)
% plot(sim_time, p_0, ...
%     sim_time, p_1, sim_time, p_2, sim_time, p_3, LineWidth=lw);
% legend("car_0", "car_1", "car_2", "car_3");
% xlabel('Time (s)')
% ylabel('Position (m)')
% grid on

v_0 = velocity.Data(:,1);
v_1 = velocity.Data(:,2);
v_2 = velocity.Data(:,3);
v_3 = velocity.Data(:,4);
v_4 = velocity.Data(:,5);
v_5 = velocity.Data(:,6);
v_6 = velocity.Data(:,7);
v_7 = velocity.Data(:,8);
v_8 = velocity.Data(:,9);
v_9 = velocity.Data(:,10);
v_10 = velocity.Data(:,11);
subplot(211)
plot(sim_time, v_0, ...
    sim_time, v_1, ...
    sim_time, v_2, ...
    sim_time, v_3, ...
    sim_time, v_4, ...
    sim_time, v_5, ...
    sim_time, v_6, ...
    sim_time, v_7, ...
    sim_time, v_8, ...
    sim_time, v_9, ...
    sim_time, v_10, ...
    LineWidth=lw); xlim([0 80]); ylim([-25 26])
% legend("car_0", "car_1", "car_2", "car_3");
ylabel('$v_i$(m/s)', 'Interpreter','latex')
grid on

fig_loc = '/home/phatdo/Dropbox/tex_doc/01_integrated_project_final_report/fig/';
fig1_name = 'sim_fixed_topo_10veh.tex';
cleanfigure('targetResolution', 100);
matlab2tikz(append(fig_loc, fig1_name), ... filename
    'width', '0.4\textwidth', ... image width
    'height', '0.3\textwidth', ... image height
     'showInfo', false);  % ... turn off information


figure(2)
subplot(212)
plot(sim_time, e_1, ...
    sim_time, e_2, ...
    sim_time, e_3, ...
    LineWidth=lw); xlim([0 80]); ylim([-5 3])
% legend("e1", "e2", "e3");
xlabel('t(s)')
ylabel('$e_i$(m)', 'Interpreter','latex')
grid on

subplot(211)
plot(sim_time, v_0, ...
    sim_time, v_1, ...
    sim_time, v_2, ...
    sim_time, v_3, ...
    LineWidth=lw); xlim([0 80]); ylim([-10 26])
% legend("car_0", "car_1", "car_2", "car_3");
ylabel('$v_i$(m/s)', 'Interpreter','latex')
grid on

% fig1_name = 'sim_fixed_topo_3veh.tex';
% cleanfigure('targetResolution', 50);
% matlab2tikz(append(fig_loc, fig1_name), ... filename
%     'width', '0.4\textwidth', ... image width
%     'height', '0.3\textwidth', ... image height
%      'showInfo', false);  % ... turn off information
