clear;
clc;

%% System definition

% Define graph
% s = [1 1 1 2 2 3];
% t = [2 3 4 3 4 4];
% G = graph(s,t);
% L = laplacian(G);
% L = eye(4)*L;
% % L = [1 0 0 -1;
% %     -1 1 0 0;
% %     0 -1 2 -1;
% %     0 -1 0 1]
% lamda_L = eig(L)
% lamda_1 = lamda_L(2);

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
Y = sdpvar(1,3,'full');
theta = sdpvar(1);

Y0 = sdpvar(1,3);
Y1 = sdpvar(1,3);
Y = Y0 + theta*Y1;

A = [theta 1 0; 0 -1-theta 1; 0 2*theta -0.3+theta]
B = [0;0;1];
C = [1 1 1];
D = [0;0;0];

% Define contraints and solve LMIs
kappa = 1;
F1 = [X*A' + A*X + lamda_1*B*Y + conj(lamda_1)*Y'*B' + 2*kappa*X <= 0];
F2 = [X*A' + A*X + lamda_3*B*Y + conj(lamda_3)*Y'*B' + 2*kappa*X <= 0];
F3 = [X*A' + A*X + lamda_4*B*Y + conj(lamda_4)*Y'*B' + 2*kappa*X <= 0];
F = [X>=0, F1, F2, F3, -2 <= theta <= 2, uncertain(theta)];

% Use an objective function
optimize(F, -trace(X));
K0 = value(Y0)*inv(value(X))
K1 = value(Y1)*inv(value(X))


%% Results
% Test the controller on simulink
sim("lpv_modeling_compact.slx");