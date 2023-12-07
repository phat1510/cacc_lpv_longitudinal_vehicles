clc
clear all
close all
%% Regression problem
x = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
A = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4+8*rand(length(t),1));
e(100:115) = 30;
y = A*x+e;
plot(t,y);

xhat = sdpvar(6,1)

%% Large-scale quadratic programs
