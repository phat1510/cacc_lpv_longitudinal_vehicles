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

tau_max = 0.12;
tau_min = 0.07;
A_min = [0 1 0; 0 0 1; 0 0 -1/tau_max];
B_min = [0 0 1/tau_max]';

A_max = [0 1 0; 0 0 1; 0 0 -1/tau_min];
B_max = [0 0 1/tau_min]';


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

%% LPV case
X = sdpvar(3,3);
Y = sdpvar(1,3,'full');
theta = sdpvar(1);
Y0 = sdpvar(1,3);
Y1 = sdpvar(1,3);
Y = Y0 + theta*Y1;
A_lpv = [0 1 0; 0 0 1; 0 0 -theta];
B_lpv = [0 0 theta]';
F1 = [X>=0];
F2 = [X*A_lpv' + A_lpv*X + lamda_2*B_lpv*Y + conj(lamda_2)*Y'*B_lpv' + 2*kappa*X<= 0];
F = [F1, F2, 8.3333 <= theta <= 14.2857, uncertain(theta)];
optimize(F,0,options);
K0 = value(Y0)*inv(value(X))
K1 = value(Y1)*inv(value(X))

%% Results

% Test on error model
% 3 vehicles
% sim("fixed_graph_sim.slx")
% 10 vehicles
sim("fixed_graph_sim_lpv.slx")


sim_time = tout;
lw = 0.5;

% Spacing error
e_1 = error.Data(:,1);
e_2 = error.Data(:,2);
e_3 = error.Data(:,3);

v_0 = velocity.Data(:,1);
v_1 = velocity.Data(:,2);
v_2 = velocity.Data(:,3);
v_3 = velocity.Data(:,4);

figure(2)
subplot(212)
plot(sim_time, e_1, ...
    sim_time, e_2, ...
    sim_time, e_3, ...
    LineWidth=lw); xlim([0 80]); %ylim([-5 3])
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
