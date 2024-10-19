function cmp_traj
clc;
clear all;
close all;

Perc = 'x';

dir = 'E:\paper2_Raoul\Sim_an_HH\LIF_IappIx\v1\';
load(strcat(dir,strcat(Perc, '.mat')))

figure(1); hold on;
wb_A_i = (18069*1e-8); % [cm^2]
plot(x_lin./wb_A_i.*(1e-12).*(1e+6), MFR_I_lin, 'r', 'LineWidth', 4);
% plot(x_lin, MFR_I_lin, 'r', 'LineWidth', 4);
grid on

dir = 'E:\paper2_Raoul\Sim_an_HH\WB_IappIx\v2\';
load(strcat(dir,strcat(Perc, '.mat')))

figure(1); hold on;
wb_A_i = (18069*1e-8); % [cm^2]
plot(x_lin./wb_A_i.*(1e-12).*(1e+6), MFR_I_lin, 'b', 'LineWidth', 4);
% plot(x_lin, MFR_I_lin, 'b', 'LineWidth', 4);
grid on

figure(1); hold on;
plot(x_lin./wb_A_i.*(1e-12).*(1e+6), x_lin./x_lin.*65, 'r', 'LineWidth', 4);
% plot(x_lin, x_lin./x_lin.*65, 'r', 'LineWidth', 4);
grid on

dir = 'E:\paper2_Raoul\Sim_an_HH\BW_IappIx\v1\';
load(strcat(dir,strcat(Perc, '.mat')))

x_low = x_lin(1:1172, 1);
MFR_I_low = MFR_I_lin(1:1172, 1);

x_high = x_lin(1173:end, 1);
MFR_I_high = MFR_I_lin(1173:end, 1);
figure(1); hold on;
wb_A_i = (18069*1e-8); % [cm^2]
plot(x_low./wb_A_i.*(1e-12).*(1e+6), MFR_I_low, 'b--', 'LineWidth', 4);
% plot(x_low, MFR_I_low, 'b--', 'LineWidth', 4);

plot(x_high./wb_A_i.*(1e-12).*(1e+6), MFR_I_high, 'b--', 'LineWidth', 4);
% plot(x_high, MFR_I_high, 'b--', 'LineWidth', 4);
% plot(MFR_I_lin, '-ko');
grid on



% set(gca,'XTick',[0 2 4 6 8 10],'XTickLabel',{'';'';'';'';'';''});
% set(gca,'YTick',[0 100 200 300],'YTickLabel',{'';'';'';''});
% 
% ylim([-10 300]);
% xlim([0 10]);
%     
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% 
% maximize_a_fig(gcf);
% savefig('cmp_freq', 'eps');

end