function my_plot
clc;clear all;close all;
format long

Perc = strcat('X');

dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEIEII\NWCA1_gIIxi_gIIGJx_0.01_20_41_IappI1280_IappE1258.9\v0\data\';
load(strcat(dir,strcat(Perc, '.mat')))

barFontSize = 25;
wb_A_i = (18069*1e-8); % [cm^2]

x_lin = unique(x_lin);
y_lin = unique(y_lin);

% x_lin_correct_unit = x_lin;
% y_lin_correct_unit = y_lin;

x_lin = 0.01*((20.0/0.01).^(x_lin/(41 - 1.0)));
x_lin_correct_unit = x_lin./wb_A_i.*(1e-9).*(1e+3);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-9).*(1e+3);

% bad_id = (30.0 > I_PowerFreq_dt);
% bad_id = (0 > PowerFreq_dt);
bad_id = (0.08 > kappa_dt(:, 2));

% E_PowerFreq_dt(bad_id) = NaN;
% E_Freq_dt(bad_id) = NaN;
I_Freq_dt(bad_id) = NaN;
% kappa_dt(bad_id, 2) = NaN;
% MFR_dt(bad_id, 1) = NaN;
% MFR_dt(bad_id, 2) = NaN;
% Freq_dt(bad_id) = NaN;
    
% tmp = reshape(E_PowerFreq_dt(:, 1), N_y, N_x);
% tmp = reshape(I_PowerFreq_dt(:, 1), N_y, N_x);
% tmp = reshape(E_Freq_dt(:, 1), N_y, N_x);
tmp = reshape(I_Freq_dt(:, 1), N_y, N_x);
% tmp = reshape(Freq_dt(:, 1), N_y, N_x);
tmp_E_kappa = reshape(kappa_dt(:, 1), N_y, N_x);
tmp_I_kappa = reshape(kappa_dt(:, 2), N_y, N_x);
tmp_E_MFR = reshape(MFR_dt(:, 1), N_y, N_x);
tmp_I_MFR = reshape(MFR_dt(:, 2), N_y, N_x);
tmp_Freq_dt = reshape(Freq_dt(:, 1), N_y, N_x);

bad_E_MFR_id = (tmp_E_MFR < 0.1);
tmp_E_kappa(bad_E_MFR_id) = NaN;
% tmp_E_MFR(bad_E_MFR_id) = NaN;
% tmp_E_MFR(bad_E_MFR_id) = NaN;

max(max(tmp'))
min(min(tmp'))

figure(1);hold on;
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 59, 'CLimEnd', 66, ...
    'is_ShowEdge', 0);
h=colorbar;

% set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
% set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('Frequency [Hz]');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on
% view(0,90);
% xlim([0.01 0.4]);
% ylim([0.001 0.01]);

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
colorbar('off');
% set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% set(h, 'fontsize', barFontSize);
% m_savefig('gII_vs_gGJII_freq_v2', 'eps');
% figure(11)
% get_colorbar(0, 59, 66, 5, 'gII_vs_gGJII_freq')

figure(2);hold on;
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_E_kappa', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 0, 'CLimEnd', 0.1, ...
    'is_ShowEdge', 0);
h=colorbar;

set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('\kappa_E');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on
% view(0,90);
% xlim([0.01 0.4]);
% ylim([0.001 0.01]);

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
colorbar('off');
set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% set(h, 'fontsize', barFontSize);
% savefig('gII_vs_gGJII_kappaE', 'eps');
% figure(21)
% get_colorbar(0, 0, 0.1, 3, 'gII_vs_gGJII_kappaE');

figure(3);hold on;
% my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_I_kappa', ...
%     'is_XLog', 1, 'is_YLog', 0, ...
%     'CLimBegin', 0, 'CLimEnd', 1.0, ...
%     'is_ShowEdge', 0);
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_I_kappa', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 0, 'CLimEnd', 0.7, ...
    'is_ShowEdge', 0);
h=colorbar;

set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('\kappa_I');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on
% view(0,90);
% xlim([0.01 0.4]);
% ylim([0.001 0.01]);

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
colorbar('off');
set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% set(h, 'fontsize', barFontSize);
% m_savefig('gII_vs_gGJII_kappaI', 'eps');
m_savefig('gII_vs_gGJII_kappaI_v2', 'eps');
% figure(31)
% get_colorbar(0, 0, 0.5, 3, 'gII_vs_gGJII_kappaI');

figure(4);hold on;
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_E_MFR', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 0, 'CLimEnd', 25, ...
    'is_ShowEdge', 0);
h=colorbar;

set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('MFR_E [spks/s]');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on
% view(0,90);
% xlim([0.01 0.4]);
% ylim([0.001 0.01]);

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
colorbar('off');
set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% set(h, 'fontsize', barFontSize);
% savefig('gII_vs_gGJII_MFRe', 'eps');
% figure(41)
% get_colorbar(0, 0, 25, 3, 'gII_vs_gGJII_MFRe');

figure(5);hold on;
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_I_MFR', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 10, 'CLimEnd', 100, ...
    'is_ShowEdge', 0);
h=colorbar;

set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('MFR_I [spks/s]');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
colorbar('off');
set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% set(h, 'fontsize', barFontSize);
% m_savefig('gII_vs_gGJII_MFRi', 'eps');
% figure(51)
% get_colorbar(0, 20, 70, 3, 'gII_vs_gGJII_MFRi');

figure(6);hold on
my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp'./tmp_I_MFR', ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 0, 'CLimEnd', 3, ...
    'is_ShowEdge', 0);
h=colorbar;

set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
% hTitle = title('MFR_I [spks/s]');
% hXLabel = xlabel('g_{I->I} [mS/cm^2]');
% hYLabel = ylabel('g_{GJ} [mS/cm^2]');
axis square
box on

make_me_pretty(gcf, ...
    gca, 15, ...
    [], 12, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
% colorbar('off');
set(gca,'XMinorGrid','Off','YMinorGrid','Off','XGrid','Off','YGrid','Off');
% m_savefig('gII_vs_gGJII_freqMFRi', 'eps');

figure(7);hold on
get_discrete_colorbar(0, 3, 59, 66, 8, 'gII_vs_gGJII_freq_v2_discrete_horizontal_bar')
end


