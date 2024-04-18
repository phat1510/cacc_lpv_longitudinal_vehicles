% Author: Minh Phat Do
% Description: Polytopic lpv cooperative control

%% Reset workspace
clear;
clc;

%% System definition

% Define directed graph
G = digraph([1 2 3 3 4], [4 1 2 4 2]);
plot(G)
Adj = adjacency(G)*eye(4);
S = [1;1;1;1];
L0 = Adj*S;
D0 = diag(L0);
L = D0-Adj
lamda_L = eig(L)
lamda_1 = lamda_L(1);
lamda_3 = lamda_L(3);
lamda_4 = lamda_L(4);

%% LMIs setup

% Define SDP variables
X = sdpvar(3,3);
Y_max = sdpvar(1,3,'full');
Y_min = sdpvar(1,3,'full');

theta_max =  2;
theta_min = -2;
A1 = [theta_min 1 0; 0 -1-theta_min 1; 0 2*theta_min -0.3+theta_min]
A2 = [theta_max 1 0; 0 -1-theta_max 1; 0 2*theta_max -0.3+theta_max]
B = [0;0;1];
C = [1 1 1];
D = [0;0;0];

% Define contraints and solve LMIs
kappa = 1;
F1 = [X*A1' + A1*X + lamda_1*B*Y_min + conj(lamda_1)*Y_min'*B' + 2*kappa*X <= 0];
F2 = [X*A2' + A2*X + lamda_1*B*Y_max + conj(lamda_1)*Y_max'*B' + 2*kappa*X <= 0];
F3 = [X*A1' + A1*X + lamda_3*B*Y_min + conj(lamda_3)*Y_min'*B' + 2*kappa*X <= 0];
F4 = [X*A2' + A2*X + lamda_3*B*Y_max + conj(lamda_3)*Y_max'*B' + 2*kappa*X <= 0];
F5 = [X*A1' + A1*X + lamda_4*B*Y_min + conj(lamda_4)*Y_min'*B' + 2*kappa*X <= 0];
F6 = [X*A2' + A2*X + lamda_4*B*Y_max + conj(lamda_4)*Y_max'*B' + 2*kappa*X <= 0];

F = [X>=0, F1, F2, F3, F4, F5, F6]; % Tr P = 1?
optimize(F);
% optimize(F, -trace(X));


%% Results
% Construct K(theta)
% K0_ref = 1e3 * [-4.1325 -1.0697 -0.1033]
% K1_ref = 1e3 * [-1.1386 -0.2887 -0.0279]
K_min = value(Y_min)*inv(value(X));
K_max = value(Y_max)*inv(value(X));
K0 = K_max + theta_max/(theta_max - theta_min)*K_min - theta_max/(theta_max - theta_min)*K_max
K1 = K_max/(theta_max - theta_min) - K_min/(theta_max - theta_min)

% K0 = K0_ref;
% K1 = K1_ref;

% Test the controller on simulink
% sim("lpv_modeling_compact.slx");