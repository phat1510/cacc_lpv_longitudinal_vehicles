% Author: Minh Phat Do
% Description: Consensus-based CACC

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
q0_init = 0;
q1_init = 0;
q2_init = 0;
q3_init = 0;

tau = 0.1;
A = [0 1 0; 0 0 1; 0 0 -1/tau];
B = [0 0 1/tau]';
C = diag([1 1 1]);
D = [0 0 0]';

error_ss = ss(A, B, C, D);
Q = diag([1 1 1]);
R = 1;
K_lqr = lqr(error_ss, Q, R)

%% LMIs setup

% Define SDP variables
X = sdpvar(3,3);
Y = sdpvar(1,3,'full');

% Define contraints and solve LMIs
kappa = 1;
F = [X>=0];
F = [F, X*A' + A*X + lamda_2*B*Y + conj(lamda_2)*Y'*B' <= 0];
options = sdpsettings();
options.verbose = 0;
optimize(F,0,options);
K_lmi = value(Y)*inv(value(X));

% K_cacc = [0.2 1.2 0] %from the author

% K_cacc = [K_lqr(1) K_lqr(2) 0];
K_cacc = K_lqr;

%% Results


% Test on error model
sim("cacc_consensus_test_full_system.slx")

sim_time = tout;

veh_0 = monitor.Data(:,1);
veh_1 = monitor.Data(:,2);
veh_2 = monitor.Data(:,3);
veh_3 = monitor.Data(:,4);

figure(1)
plot(sim_time, veh_0, sim_time, veh_1, ...
    sim_time, veh_2, sim_time, veh_3, LineWidth=2);
legend("car_0", "car_1", "car_2", "car_3");
ylabel('Error(m)')
grid on
