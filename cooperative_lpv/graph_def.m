%% Reset workspace
clear
clc
close

%%
% s = [1 1 1 2 2 3];
% t = [2 3 4 3 4 4];
% G = graph(s,t);
G = digraph([1 2 3 3 4], [4 1 2 4 2])
plot(G)
A = adjacency(G)*eye(4);
B = [1;1;1;1];
L0 = A*B;
D = diag(L0);
L = D-A
E = eig(L)
