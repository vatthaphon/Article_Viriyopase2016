function my_plot_copy
clc;clear all;
close all;
format long

nw_A_e = (21590*1e-8); % [cm^2]
wb_A_i = (18069*1e-8); % [cm^2]
factor1 = (1e-12);
factor2 = (1e+6);

% nw_A_e = 1; % [cm^2]
% wb_A_i = 1; % [cm^2]
% factor1 = 1;
% factor2 = 1;

Perc = 'X';
Patch_FaceAlpha = 0.4;
PINGING_FaceAlpha = 0.5 + 0.2;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.4;

Ii_min = 0;
Ii_max = 10000;

Fixed_Ie = 17;
Fixed_Ii = 7.5;

begin_Val13 = 10;
end_Val13 = 10000;
len_Val13 = 41;

label_fz = 20;
gca_fz = 20;

%% PING + ING
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEEIIEII\NWCA1_IappExi_IappIx_10_10000_41_sigmaWNE60_sigmaWNI0.5\v0\data\';
load(strcat(dir,strcat(Perc, '.mat')))

kappa_dt_for_PING = kappa_dt;
PowerFreq_dt_for_PING = PowerFreq_dt;
Freq_dt_for_PING = Freq_dt;
E_Freq_dt_for_PING = E_Freq_dt;
I_Freq_dt_for_PING = I_Freq_dt;
MFR_dt_for_PING = MFR_dt;

N_x_for_PING = size(unique(x_lin), 1);
N_y_for_PING = size(unique(y_lin), 1);

good_y_lin = (y_lin >= Ii_min) & (y_lin <= Ii_max) ;
x_lin = x_lin(good_y_lin);
y_lin = y_lin(good_y_lin);
E_Freq_dt = E_Freq_dt(good_y_lin, 1);
I_Freq_dt = I_Freq_dt(good_y_lin, 1);
E_PowerFreq_dt = E_PowerFreq_dt(good_y_lin, 1);
I_PowerFreq_dt = I_PowerFreq_dt(good_y_lin, 1);
E_MFR_dt = MFR_dt(good_y_lin, 1);
I_MFR_dt = MFR_dt(good_y_lin, 2);
E_kappa_dt = kappa_dt(good_y_lin, 1);
I_kappa_dt = kappa_dt(good_y_lin, 2);

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

% Create a 3D matrix
x_lin = unique(x_lin);
y_lin = unique(y_lin);

x_lin = begin_Val13*((end_Val13/begin_Val13).^(x_lin/(len_Val13 - 1.0)));

% x_lin_correct_unit = x_lin;
% y_lin_correct_unit = y_lin;

x_lin_correct_unit = x_lin./nw_A_e.*factor1.*factor2;
y_lin_correct_unit = y_lin./wb_A_i.*factor1.*factor2;

x_lin_correct_unit_tmp = x_lin_correct_unit;

x_lin_correct_unit = log10(x_lin_correct_unit);

N_x = size(x_lin_correct_unit, 1);
N_y = size(y_lin_correct_unit, 1);

b_rhythms = (E_PowerFreq_dt < 0.01) | (isnan(E_PowerFreq_dt)) | (isnan(MFR_dt(:, 1)))  | (MFR_dt(:, 1) <= 0.15) ;
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;

figure(4);hold on;
tmp = reshape(E_Freq_dt, N_y, N_x);
z_const2 = tmp(41, :);
max_y = max(y_lin_correct_unit);
plot3(x_lin_correct_unit, ones(size(x_lin_correct_unit), 1)*max_y, z_const2, 'LineWidth', 20, 'Color', [0 0 0]);
hSurface1 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
vary_Ii_PINGING_freq = tmp(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_PINGING_freq = tmp(id, :);

% figure(41);hold on;
% z_const1 = tmp(1, :);
% z_const2 = tmp(11, :);
% z_const3 = tmp(:, 1);
% z_const4 = tmp(:, 11);
% plot(x_lin_correct_unit, z_const2);

figure(5);hold on;
tmp = reshape(E_MFR_dt, N_y, N_x);
Fig5_PINGING_x_lin_correct_unit = x_lin_correct_unit;
Fig5_PINGING_y_lin_correct_unit = y_lin_correct_unit;
Fig5_PINGING_tmp = tmp;

figure(6);hold on;
tmp = reshape(I_MFR_dt, N_y, N_x);
hSurface3 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

% id = sum(x_lin_correct_unit <= Fixed_Ie);
% vary_Ii_PINGING_freq = tmp(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_PINGING_freq = tmp(id, :);

% figure(5);hold on;
% tmp = reshape(MFR_dt(:, 1), N_y, N_x);
% % tmp = reshape(E_PowerFreq_dt, N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
% 
% figure(6);hold on;
% tmp = reshape(MFR_dt(:, 2), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
% 
% figure(7);hold on;
% tmp = reshape(kappa_dt(:, 1), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
% 
% figure(8);hold on;
% tmp = reshape(kappa_dt(:, 2), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

%% PING
Perc1 = 'X';
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEEIIEII\NWCA1_IappExi_10_10000_41_IappI0_sigmaWNE60_sigmaWNI0.5\v0\data\';
load(strcat(dir, strcat(Perc1, '.mat')))

running_i = 1;
for i = 1:N_x_for_PING
    E_kappa_dtMAT(i,1) = kappa_dt(running_i, 1);
    I_kappa_dtMAT(i,1) = kappa_dt(running_i, 2);
    PowerFreq_dtMAT(i,1) = PowerFreq_dt(running_i, 1);
    Freq_dtMAT(i,1) = Freq_dt(running_i, 1);
    E_Freq_dtMAT(i,1) = E_Freq_dt(running_i, 1);
    I_Freq_dtMAT(i,1) = I_Freq_dt(running_i, 1);
    E_MFR_dtMAT(i,1) = MFR_dt(running_i, 1);
    I_MFR_dtMAT(i,1) = MFR_dt(running_i, 2);
    
    running_i = running_i + 1;
end

for i = 1:1:N_x_for_PING
    for j = 1:1:size(y_lin_correct_unit, 1)
        freq(i, j) = Freq_dtMAT(i, 1);
        e_freq(i, j) = E_Freq_dtMAT(i, 1);
        e_mfr(i, j) = E_MFR_dtMAT(i, 1);
        i_mfr(i, j) = I_MFR_dtMAT(i, 1);
        e_kappa(i, j) = E_kappa_dtMAT(i, 1);
        i_kappa(i, j) = I_kappa_dtMAT(i, 1);        
    end
end

figure(4);hold on;
tmp = e_freq';
hSurface4 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
vary_Ii_PING_freq = tmp(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_PING_freq = tmp(id, :);

figure(5);hold on;
tmp = e_mfr';
Fig5_PING_x_lin_correct_unit = x_lin_correct_unit;
Fig5_PING_y_lin_correct_unit = y_lin_correct_unit;
Fig5_PING_tmp = tmp;

figure(6);hold on;
tmp = i_mfr';
hSurface6 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

% id = sum(x_lin_correct_unit <= Fixed_Ie);
% vary_Ii_PING_freq = tmp(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_PING_freq = tmp(id, :);

% figure(5);hold on;
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, e_mfr');
% set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
% 
% figure(6);hold on;
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, i_mfr');
% set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
% 
% figure(7);hold on;
% tmp = reshape(e_kappa, N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
% 
% figure(8);hold on;
% tmp = reshape(i_kappa, N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);

%% ING
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEIEII\NWCA1_IappExi_IappIx_10_10000_41_sigmaWNE60_sigmaWNI0.5\v0\data\';
load(strcat(dir,strcat(Perc, '.mat')))

% kappa_dt_for_PING = kappa_dt;
% PowerFreq_dt_for_PING = PowerFreq_dt;
% Freq_dt_for_PING = Freq_dt;
% E_Freq_dt_for_PING = E_Freq_dt;
% I_Freq_dt_for_PING = I_Freq_dt;
% MFR_dt_for_PING = MFR_dt;
% 
% N_x_for_PING = size(unique(x_lin), 1);
% N_y_for_PING = size(unique(y_lin), 1);

good_y_lin = (y_lin >= Ii_min) & (y_lin <= Ii_max);
x_lin = x_lin(good_y_lin);
y_lin = y_lin(good_y_lin);
E_Freq_dt = E_Freq_dt(good_y_lin, 1);
I_Freq_dt = I_Freq_dt(good_y_lin, 1);
E_PowerFreq_dt = E_PowerFreq_dt(good_y_lin, 1);
I_PowerFreq_dt = I_PowerFreq_dt(good_y_lin, 1);
E_MFR_dt = MFR_dt(good_y_lin, 1);
I_MFR_dt = MFR_dt(good_y_lin, 2);
E_kappa_dt = kappa_dt(good_y_lin, 1);
I_kappa_dt = kappa_dt(good_y_lin, 2);

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

% Create a 3D matrix
% x_lin = unique(x_lin);
% y_lin = unique(y_lin);
% 
% % x_lin_correct_unit = x_lin;
% % y_lin_correct_unit = y_lin;
% 
% x_lin = begin_Val13*((end_Val13/begin_Val13).^(x_lin/(len_Val13 - 1.0)));
% 
% x_lin_correct_unit = x_lin./nw_A_e.*factor1.*factor2;
% y_lin_correct_unit = y_lin./wb_A_i.*factor1.*factor2;

N_x = size(x_lin_correct_unit, 1);
N_y = size(y_lin_correct_unit, 1);

b_rhythms = (E_PowerFreq_dt < 0.003) | (isnan(E_PowerFreq_dt));
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;
I_Freq_dt(b_I_rhythms, 1) = NaN;

figure(4);hold on;
tmp = reshape(I_Freq_dt, N_y, N_x);
hSurface7 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
min_ING_freq = min(min(tmp));
max_ING_freq = max(max(tmp));
min_ING_x = min(min(x_lin_correct_unit));
max_ING_x = max(max(x_lin_correct_unit));
min_ING_y = min(min(y_lin_correct_unit));
max_ING_y = max(max(y_lin_correct_unit));

id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
vary_Ii_ING_freq = tmp(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_ING_freq = tmp(id, :);

figure(5);hold on;
tmp = reshape(MFR_dt(:, 1), N_y, N_x);
Fig5_ING_x_lin_correct_unit = x_lin_correct_unit;
Fig5_ING_y_lin_correct_unit = y_lin_correct_unit;
Fig5_ING_tmp = tmp;

figure(6);hold on;
tmp = reshape(MFR_dt(:, 2), N_y, N_x);
hSurface9 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

% id = sum(x_lin_correct_unit <= Fixed_Ie);
% vary_Ii_ING_freq = tmp(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_ING_freq = tmp(id, :);

% figure(5);hold on;
% tmp = reshape(MFR_dt(:, 1), N_y, N_x);
% % tmp = reshape(I_PowerFreq_dt, N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
% 
% figure(6);hold on;
% tmp = reshape(MFR_dt(:, 2), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
% 
% figure(7);hold on;
% tmp = reshape(kappa_dt(:, 1), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
% 
% figure(8);hold on;
% tmp = reshape(kappa_dt(:, 2), N_y, N_x);
% hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

min_log = floor(min(x_lin_correct_unit));
max_log = ceil(max(x_lin_correct_unit));

min_log_val = str2num(strcat('10^', num2str(min_log)));
max_log_val = str2num(strcat('10^', num2str(max_log)));

id = 1;
x = min_log_val;
XTickArray = [];
while x < max_log_val
    x = min_log_val*id;
    XTickArray = [XTickArray x];
    id = id + 1;
    
    if (id == 10)
        id = 1;
        min_log_val = min_log_val*10;
    end
end

[rows, cols] = size(XTickArray);
XTickLabelArray = cell(rows, cols);

id = 1;
j = 1;
for i = XTickArray
%     if (id == 1)
%         XTickLabelArray{j} = strcat('10^', num2str(log10(i)));
%     else
        XTickLabelArray{j} = '';
%     end
    id = id + 1;
    j = j + 1;
    
    if (id == 10)
        id = 1;
    end
end

XTickArray = log10(XTickArray);

figure(4);hold on;
box on
% xlabel('I_{0,E} [\muA/cm^2]');
% ylabel('I_{0,I} [\muA/cm^2]');
% zlabel('Freq. [Hz]');
axis tight
axis square
view([-17.5, 4])
% set(gca, 'Xscale', 'log');

make_me_pretty(gcf, ...
    gca, gca_fz, ...
    [], 12, ...
    get(gca,'xlabel'), label_fz, ...
    get(gca,'ylabel'), label_fz, ...
    get(gca,'zlabel'), label_fz, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',XTickArray,'XTickLabel',XTickLabelArray);
set(gca, 'XMinorTick'  , 'off');

% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.6 0.8 1 1.2 1.4],'YTickLabel',{'';'';'';'';''});
% set(gca,'ZTick',[30 35 40 45],'ZTickLabel',{'';'';'';''});

% maximize_a_fig(gcf);
% savefig('PINGING_IeIi_Efreq', 'eps');

figure(5);hold on;
Fig5_PINGING_tmp(abs(Fig5_ING_tmp - Fig5_PINGING_tmp) <= 0.2) = NaN;

hSurface2 = surf(Fig5_PINGING_x_lin_correct_unit, Fig5_PINGING_y_lin_correct_unit, Fig5_PINGING_tmp);
hSurface5 = surf(Fig5_PING_x_lin_correct_unit, Fig5_PING_y_lin_correct_unit, Fig5_PING_tmp + 1);
hSurface8 = surf(Fig5_ING_x_lin_correct_unit, Fig5_ING_y_lin_correct_unit, Fig5_ING_tmp);

box on
% xlabel('I_{0,E} [\muA/cm^2]');
% ylabel('I_{0,I} [\muA/cm^2]');
% zlabel('MFR_E [spks/s]');
axis tight
axis square
view([-17.5, 4])
% set(gca, 'Xscale', 'log');

make_me_pretty(gcf, ...
    gca, gca_fz, ...
    [], 12, ...
    get(gca,'xlabel'), label_fz, ...
    get(gca,'ylabel'), label_fz, ...
    get(gca,'zlabel'), label_fz, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',XTickArray,'XTickLabel',XTickLabelArray);
set(gca, 'XMinorTick'  , 'off');

% maximize_a_fig(gcf);
% savefig('PINGING_IeIi_MFRe', 'eps');

figure(6);hold on;
box on
% xlabel('I_{0,E} [\muA/cm^2]');
% ylabel('I_{0,I} [\muA/cm^2]');
% zlabel('MFR_I [spks/s]');
axis tight
axis square
view([-17.5, 4])
% set(gca, 'Xscale', 'log');

make_me_pretty(gcf, ...
    gca, gca_fz, ...
    [], 12, ...
    get(gca,'xlabel'), label_fz, ...
    get(gca,'ylabel'), label_fz, ...
    get(gca,'zlabel'), label_fz, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',XTickArray,'XTickLabel',XTickLabelArray);
set(gca, 'XMinorTick'  , 'off');

set(hSurface1,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
set(hSurface2,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
set(hSurface3,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

set(hSurface4,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
set(hSurface5,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
set(hSurface6,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);

set(hSurface7,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
set(hSurface8,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
set(hSurface9,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

figure(4)
set(gca, 'YTick', [6.6 6.8 7 7.2 7.4 7.6], 'YTickLabel',{'';'';'';'';'';''});
set(gca, 'ZTick', [50 55 60 65 70 75], 'ZTickLabel',{'';'';'';'';'';''});
maximize_a_fig(gcf);
% savefig('PINGING_IeIi_Efreq', 'eps');

figure(5)
set(gca, 'YTick', [6.6 6.8 7 7.2 7.4 7.6], 'YTickLabel',{'';'';'';'';'';''});
set(gca, 'ZTick', [20 40 60 80 100 120], 'ZTickLabel',{'';'';'';'';'';''});
maximize_a_fig(gcf);
% savefig('PINGING_IeIi_MFRe', 'eps');

figure(6)
set(gca, 'YTick', [6.6 6.8 7 7.2 7.4 7.6], 'YTickLabel',{'';'';'';'';'';''});
set(gca, 'ZTick', [10 20 30 40 50 60 70], 'ZTickLabel',{'';'';'';'';'';'';''});
maximize_a_fig(gcf);
% savefig('PINGING_IeIi_MFRi', 'eps');

% maximize_a_fig(gcf);
% savefig('PINGING_IeIi_MFRi', 'eps');

% figure(5);hold on;
% % xlabel('I_{E}');
% % ylabel('I_{I}');
% axis tight
% axis square
% % colorbar
% view([-37.500000000000000, 30])
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
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.6 0.8 1 1.2],'YTickLabel',{'';'';'';''});
% set(gca,'ZTick',[10 20 30 40 50],'ZTickLabel',{'';'';'';'';''});
% 
% % savefig('PINGING_MFRe', 'eps');
% 
% figure(6);hold on;
% % xlabel('I_{E}');
% % ylabel('I_{I}');
% axis tight
% axis square
% % colorbar
% view([-37.500000000000000, 30])
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
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.6 0.8 1 1.2],'YTickLabel',{'';'';'';''});
% set(gca,'ZTick',[20 25 30 35 40],'ZTickLabel',{'';'';'';'';''});
% 
% % savefig('PINGING_MFRi', 'eps');
% 
% % figure(7);hold on;
% % xlabel('I_{E}');
% % ylabel('I_{I}');
% % axis tight
% % axis square
% % colorbar
% % view([-37.500000000000000, 30])
% % 
% % figure(8);hold on;
% % xlabel('I_{E}');
% % ylabel('I_{I}');
% % axis tight
% % axis square
% % colorbar
% % view([-37.500000000000000, 30])

%% Create the cross section
% % X1 = [Fixed_Ie min_ING_y min_ING_freq];
% % X2 = [Fixed_Ie max_ING_y min_ING_freq];
% % X3 = [Fixed_Ie max_ING_y max_ING_freq];
% % X4 = [Fixed_Ie min_ING_y max_ING_freq];
% % makePatches(X1, X2, X3, X4, Patch_FaceAlpha);
% % 
% % X1 = [min_ING_x Fixed_Ii min_ING_freq];
% % X2 = [max_ING_x Fixed_Ii min_ING_freq];
% % X3 = [max_ING_x Fixed_Ii max_ING_freq];
% % X4 = [min_ING_x Fixed_Ii max_ING_freq];
% % makePatches(X1, X2, X3, X4, Patch_FaceAlpha);
% 
% figure(100);hold on;
% plot(y_lin_correct_unit, vary_Ii_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
% plot(y_lin_correct_unit, vary_Ii_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
% plot(y_lin_correct_unit, vary_Ii_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
% 
% axis tight
% axis square
% box on
% % ylim([min_ING_freq max_ING_freq]);
% ylim([30 46]);
% 
% % set(gca,'XTick',[0.6 0.8 1 1.2 1.4],'XTickLabel',{'';'';'';'';''});
% % set(gca,'YTick',[30 35 40 45],'YTickLabel',{'';'';'';''});
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
% % savefig('vary_Ii_freq_SigmaI0_5_for_paper2_Fixed_Ie_2', 'eps');
% 
% figure(101);hold on;
% plot(x_lin_correct_unit, vary_Ie_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
% plot(x_lin_correct_unit, vary_Ie_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
% plot(x_lin_correct_unit, vary_Ie_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
% 
% axis tight
% axis square
% box on
% % ylim([min_ING_freq max_ING_freq]);
% ylim([30 46]);
% 
% % set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% % set(gca,'YTick',[30 35 40 45],'YTickLabel',{'';'';'';''});
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
% % savefig('vary_Ie_freq_SigmaI0_5_for_paper2_Fixed_Ii_0_85', 'eps');

close all

% figure(100);hold on;
% plot(x_lin_correct_unit, vary_Ie_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
% plot(x_lin_correct_unit, vary_Ie_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
% plot(x_lin_correct_unit, vary_Ie_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
% 
% axis tight
% axis square
% set(gca,'XTick',XTickArray,'XTickLabel',XTickLabelArray);
% set(gca, 'XMinorTick'  , 'off');
% 
% set(gca,'YTick', [45 55 65 75],'YTickLabel', {''; ''; ''; ''});
% 
% set(gcf, 'color', 'white');
% set(gca, ...
%   'box'         , 'on'      , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XGrid'       , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'ZGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'ZColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1         );
% 
% box on
% ylim([42 80]);
% 
% maximize_a_fig(gcf);
% savefig('vary_Ie_freq', 'eps');

figure(101);hold on;
plot(y_lin_correct_unit, vary_Ii_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
plot(y_lin_correct_unit, vary_Ii_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');

axis tight
axis square
set(gca, 'XTick', [6.6 6.8 7.0 7.2 7.4 7.6], 'XTickLabel', [''; ''; ''; ''; ''; '';]);
set(gca, 'XMinorTick'  , 'off');

set(gca,'YTick', [56 60 64 68],'YTickLabel', {''; ''; ''; ''});

set(gcf, 'color', 'white');
set(gca, ...
  'box'         , 'on'      , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'ZGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'ZColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

box on
ylim([55 69]);

maximize_a_fig(gcf);
% savefig('vary_Ii_freq', 'eps');

end

