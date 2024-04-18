% Author: Minh Phat Do
% Description: Consensus-based CACC with 
% dynamic communication topologies

%% Reset workspace
clear; clc; close all;

%% System definition

% Define the plattoon configuration
% G = digraph([1 2 3], [2 3 4]); look-forward
% G = digraph([1 2 3 1 2 1], [2 3 4 3 4 4]);
G = digraph([1 3 2 4 3], [2 2 3 3 4]);
Adj = adjacency(G)*eye(4);
Adj = Adj'
S = [1;1;1;1];
L0 = Adj*S;
D0 = diag(L0);
L = D0-Adj;
P = diag([0 0 0 0]);
L_hat = L + P
lamda_L_hat = eig(L_hat)
lamda_1 = lamda_L_hat(1);
lamda_2 = lamda_L_hat(2);
lamda_3 = lamda_L_hat(3);

L1 = [0 0 0 0;
    -1 1 0 0;
    0 -1 1 0;
    0 0 -1 1];
L2 = [0 0 0 0;
    -1 2 -1 0;
    0 -1 2 -1;
    0 0 -1 1];

%% System definition
lv = 4.46;
h = 1;
rv = 0;
q0_init = 0;
q1_init = 0;
q2_init = 0;
q3_init = 0;

tau = 0.1;
A = [0 1 0; 0 0 1; 0 0 -1/tau];
B = [0 0 1/tau]';
C = diag([1 1 1]);
D = [0 0 0]';

% error_ss = ss(A, B, C, D);
% Q = diag([1 1 1]);
% R = 1;
% K_lqr = lqr(error_ss, Q, R);

%% LMIs setup

% Define SDP variables
X = sdpvar(3,3);
Y = sdpvar(1,3,'full');

% Define contraints and solve LMIs
kappa = 1;
F1 = [X*A' + A*X + lamda_1*B*Y + conj(lamda_1)*Y'*B' + 2*kappa*X<= 0];
F2 = [X*A' + A*X + lamda_2*B*Y + conj(lamda_2)*Y'*B' + 2*kappa*X<= 0];
F3 = [X*A' + A*X + lamda_3*B*Y + conj(lamda_3)*Y'*B' + 2*kappa*X<= 0];
F = [X>=0, F1, F2, F3];
options = sdpsettings();
options.verbose = 0;
% optimize(F,0,options);
optimize(F,0,options);
K_lmi = value(Y)*inv(value(X));

% K_cacc = [0.2 1.2 0] %from the author

% K_cacc = [K_lqr(1) K_lqr(2) 0];
K_cacc = -K_lmi
K_cacc1 = [2.3231, 2.3783, -0.3394]
K_cacc2 = [ 19.3466, 17.8329, 1.8423]
% K_cacc1 = [0 0 0];
% K_cacc2 = [0 0 0];

%% Results


% Test on error model
sim("cacc_consensus_dyn_graph_sim.slx")

sim_time = tout;

e_1 = monitor.Data(:,1);
e_2 = monitor.Data(:,2);
e_3 = monitor.Data(:,3);

lw = 1;

figure
subplot(221)
plot(sim_time, e_1, ...
    sim_time, e_2, sim_time, e_3, LineWidth=lw);
legend("1-0", "2-1", "3-2");
xlabel('Time (s)')
ylabel('Error (m)')
grid on

p_0 = monitor.Data(:,4);
p_1 = monitor.Data(:,5);
p_2 = monitor.Data(:,6);
p_3 = monitor.Data(:,7);

subplot(222)
plot(sim_time, p_0, ...
    sim_time, p_1, sim_time, p_2, sim_time, p_3, LineWidth=lw);
legend("car_0", "car_1", "car_2", "car_3");
legend('Location','northwest');
xlabel('Time (s)')
ylabel('Position (m)')
grid on

v_0 = monitor.Data(:,8);
v_1 = monitor.Data(:,9);
v_2 = monitor.Data(:,10);
v_3 = monitor.Data(:,11);
subplot(223)
plot(sim_time, v_0, ...
    sim_time, v_1, sim_time, v_2, sim_time, v_3, LineWidth=lw);
legend("car_0", "car_1", "car_2", "car_3");
legend('Location','northwest');
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

subplot(224)
plot(G)

% fig_loc = '/home/phatdo/Dropbox/tex_doc/mobile_robotics_uav_report/figures/';
% fig1_name = '9-3.fig2.tex';
% cleanfigure('targetResolution', 100);
% matlab2tikz(append(fig_loc, fig1_name), ... filename
%     'width', '0.6\textwidth', ... image width
%     'height', '0.15\textwidth', ... image height
%      'showInfo', false);  % ... turn off information
