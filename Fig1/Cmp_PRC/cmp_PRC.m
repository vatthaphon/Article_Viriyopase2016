function cmp_PRC
clc;
clear all;
close all;

N_step = 10;
Perc = 'PRC_Borger_Walker';

PRC = [];
phase = [];
dir = 'C:\paper2_Raoul\Sim_an_HH\Cmp_PRC\';
load(strcat(dir,strcat(Perc, '.mat')))

cols = size(PRC, 2);

max_PRC = max(max(PRC(1, :)));
figure(1); hold on;
good_phase = phase([1:N_step:cols]);
good_PRC = PRC(1:N_step:cols)/max_PRC;
plot([good_phase phase(end)], [good_PRC PRC(end)/max_PRC], 'r--', 'LineWidth', 4);

Perc = 'PRC_Wang_Buzsaki';

dir = 'C:\paper2_Raoul\Sim_an_HH\Cmp_PRC\';
load(strcat(dir,strcat(Perc, '.mat')))

max_PRC = max(max(PRC(1, :)));
figure(1); hold on;
good_phase = phase([1:N_step:cols]);
good_PRC = PRC(1:N_step:cols)/max_PRC;
plot([good_phase phase(end)], [good_PRC PRC(end)/max_PRC], 'r-', 'LineWidth', 4);

set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],'XTickLabel',{'';'';'';'';'';''});
set(gca,'YTick',[0 0.4 0.8 1.2],'YTickLabel',{'';'';'';''});

ylim([-0.2 1.2]);
xlim([0 1]);
    
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
% savefig('cmp_PRC_test', 'eps');

end