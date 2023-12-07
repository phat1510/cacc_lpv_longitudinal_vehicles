clc;
% clear;
close all;
B = [0;0;1];
C = [1 1 1];
D = [0 ;0; 0];

plot_on = false;

if (plot_on)
    figure(1)
    subplot(2,1,1)
    t = out.t.Time;
    theta = out.theta.Data;
    plot(t,theta);
%     xlabel('test'+'$\theta(t)$','interpreter','latex','FontSize',12);
    ylabel('$\theta(t)$','interpreter','latex','FontSize',12);
    grid on

    subplot(2,1,2)
    x1 = out.first_state.Data(:,1);
    x2 = out.first_state.Data(:,2);
    x3 = out.first_state.Data(:,3);
    x4 = out.first_state.Data(:,4);
    plot(t, x1, t, x2, t, x3, t, x4);
    ylabel('$x_k^{(1)}(t)$','interpreter','latex','FontSize',12);
    legend("Agent 1","Agent 2","Agent 3","Agent 4",'interpreter','latex','FontSize',12)
    grid on
    matlab2tikz('homo.tex','width', '12cm', 'showInfo', false);
end