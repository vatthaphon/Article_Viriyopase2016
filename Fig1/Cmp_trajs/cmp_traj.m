function cmp_traj
clc;
clear all;
close all;

figure(1); hold on;
load('E:\paper2_Raoul\Sim_an_HH\Traj_WB_I_1_1.mat');
plot(tc(:, 1), lc(:, 1), 'b', 'LineWidth', 4);

load('E:\paper2_Raoul\Sim_an_HH\Traj_BW_I_7_1.mat');
plot(tc(:, 1), lc(:, 1), 'b--', 'LineWidth', 4);

% load('E:\paper2_Raoul\Sim_an_HH\Traj_Nowacki_I_2.mat');
% plot(tc(:, 1), lc(:, 1), 'r', 'LineWidth', 4);

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[960 970 980 990],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[-80 -40 0 40],'YTickLabel',{'';'';'';''});

% legend('', '')

ylim([-100 60]);
xlim([960 990]);
grid on;

maximize_a_fig(gcf);

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

savefig('cmp_traj', 'eps');

end