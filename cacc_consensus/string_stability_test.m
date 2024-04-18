clc
clear
close all

h = 1;
tau = 0.1;

K_cacc1 = [2.3231, 2.3783, -0.3394];

s = tf('s');
G = 1/(s^2*(tau*s+1));
H = h*s+1;
K = K_cacc1(1) + K_cacc1(2)*s;

ST = 1/H * (K*G + 1)/(1+K*G);
[sv,wout] = sigma(ST);
close

% -----------------------------
figure
semilogx(wout,20*log10(sv))
xlabel('(rad/s)'), ylabel('(dB)')
grid on
% -----------------------------
fig_loc = '/home/phatdo/Dropbox/tex_doc/integrated_project_3rd_presentation/tikz/';
fig1_name = 'string_stability.tex';
cleanfigure('targetResolution', 100);
matlab2tikz(append(fig_loc, fig1_name), ... filename
    'width', '0.5\textwidth', ... image width
    'height', '0.2\textwidth', ... image height
     'showInfo', false);  % ... turn off information
% -----------------------------