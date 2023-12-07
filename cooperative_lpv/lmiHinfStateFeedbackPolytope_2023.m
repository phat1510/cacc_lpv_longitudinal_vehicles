% Author: Olivier Sename
% Nov 2023
% 
% Description
% Function that compute the state feedback Hinf controller solving the LMI problem. 
% To use this function one have to create the 
% generalized plant P such that:
%  Xdot    A  Bw  Bu    X
%   Z   =  Cz Dzw Dzu   W
% 
% Input
%  listP  : list of plants
%  ncon   : number of control signals
%  nstate : number of state variables
%  sp     : the controller is striclty proper (1=>yes, 0=>no)
%  percentage : percentage added to the gamma optimal to conditionate the
%  controller

% Output
%  listK  : list of controller (state space)
%  listCL : list of closed loop (state space)
%  gopt   : optimal gamma
%
% [listK,listCL,gopt] = lmiHinfPolytope(listP,nmeas,ncon,sp,percentage,solver)

function [K,gopt] = lmiHinfStateFeedbackPolytope_2023(listP,nstate,ncon,sp,percentage,solver)
%%% Size of the generalized plant
sizeX = size(listP{1}.a,1);
sizeZ = size(listP{1},1); % number of controlled output 
sizeW = size(listP{1},2)-ncon;  % number of input - number of control
sizeU = ncon;

%%% To tackle with some strict inequalities problems not well implemented 
%%% in Yalmip
epsi = 1e-6;

%%% Shows vector lenght
input  = [sizeX sizeW sizeU];
output = [sizeX sizeZ];

%%% Cut the system 
for i = 1:size(listP,2)
    A{i}   = listP{i}.a(1:sizeX,1:sizeX);
    Bw{i}  = listP{i}.b(1:sizeX,1:sizeW);
    Bu{i}  = listP{i}.b(1:sizeX,sizeW+1:sizeW+sizeU);
    Cz{i}  = listP{i}.c(1:sizeZ,1:sizeX);
    Dzw{i} = listP{i}.d(1:sizeZ,1:sizeW);
    Dzu{i} = listP{i}.d(1:sizeZ,sizeW+1:sizeW+sizeU);

end;

% %%% Create variables matrix
for i = 1:size(listP,2)
    Y{i}     = sdpvar(sizeU,sizeX); 
end;
X     = sdpvar(sizeX,sizeX,'symmetric'); 
gamma = sdpvar(1,1,'full');%1;%

%%% LMIs definition of the Hinf problem
F=[X>=epsi];
%F  = set(H0>eps,'X and Y constraint');
for i = 1:size(listP,2)
    M11{i}=A{i}*X+X*A{i}'+Bu{i}*Y{i}+(Bu{i}*Y{i})'; 
    M21{i}=Bw{i}';
    M31{i}=Cz{i}*X+Dzu{i}*Y{i};
    M32{i}=Dzw{i};
    H1{i}  = [ M11{i}        M21{i}'               M31{i}';
               M21{i}      -eye(sizeW)           M32{i}';
               M31{i}        M32{i}          -gamma*eye(sizeZ)]; 
    name = ['Hinf constraint ' int2str(i)];
    F = [F,H1{i}<=-epsi];
end;
%F

%%% Find feasible solution minimizing gamma
%solver='sedumi';
ops      = sdpsettings('solver',solver);
solution = optimize(F,gamma,ops);

gopt = sqrt(double(gamma))*(1+percentage/100);

% % 
% % %%% Extract numerical solution
for i = 1:size(listP,2)
Y{i}=double(Y{i});
end;
X = double(X);

for i = 1:size(listP,2)
K{i}= Y{i}*inv(X);
end;

%%% Some verifications
VPx = eig(X)
for i=1:size(X,1)
    if (VPx(i) <= 0)
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('Error, Lyapunov function are non positive definite')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end;
end;
