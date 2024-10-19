function my_plot
clc;clear all;
close all;
format long

Perc = 'X';
Patch_FaceAlpha = 0.4;
PINGING_FaceAlpha = 0.5 + 0.0;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.0;
Ii_min = 90;
Ii_max = 270;

Fixed_Ie = 2.0;
Fixed_Ii = 0.85;
% Fixed_Ie = 0 + 1000;

Fixed_Ie5 = 2.0;
Fixed_Ie4 = 1.5;
Fixed_Ie3 = 1.0;
Fixed_Ie2 = 0.5;
Fixed_Ie1 = 0.0;

% Fixed_Ie = 432.0;
% Fixed_Ii = 0.85 + 1000;

nw_A_e = (21590*1e-8); % [cm^2]
wb_A_i = (18069*1e-8); % [cm^2]

%% PING + ING
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0.3_pII0.3_pIIGJ0.004\v6\';
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

% x_lin_correct_unit = x_lin;
% y_lin_correct_unit = y_lin;

x_lin_correct_unit = x_lin./nw_A_e.*(1e-12).*(1e+6);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-12).*(1e+6);

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
tmp = reshape(E_Freq_dt, N_y, N_x); % I use this by default.
% tmp = reshape(I_Freq_dt, N_y, N_x); % This is not what I use to get the figure in the paper.
tmp_E_MFR_dt = reshape(E_MFR_dt, N_y, N_x);
tmp_I_MFR_dt = reshape(I_MFR_dt, N_y, N_x);
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp + 0.1);
set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

id = sum(x_lin_correct_unit <= Fixed_Ie);
vary_Ii_PINGING_freq = tmp(:, id);
vary_Ii_PINGING_E_MFR = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR = tmp_I_MFR_dt(:, id);

id = sum(x_lin_correct_unit <= Fixed_Ie1);
vary_Ii_PINGING_freq1 = tmp(:, id);
vary_Ii_PINGING_E_MFR1 = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR1 = tmp_I_MFR_dt(:, id);

id = sum(x_lin_correct_unit <= Fixed_Ie2);
vary_Ii_PINGING_freq2 = tmp(:, id);
vary_Ii_PINGING_E_MFR2 = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR2 = tmp_I_MFR_dt(:, id);

id = sum(x_lin_correct_unit <= Fixed_Ie3);
vary_Ii_PINGING_freq3 = tmp(:, id);
vary_Ii_PINGING_E_MFR3 = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR3 = tmp_I_MFR_dt(:, id);

id = sum(x_lin_correct_unit <= Fixed_Ie4);
vary_Ii_PINGING_freq4 = tmp(:, id);
vary_Ii_PINGING_E_MFR4 = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR4 = tmp_I_MFR_dt(:, id);

id = sum(x_lin_correct_unit <= Fixed_Ie5);
vary_Ii_PINGING_freq5 = tmp(:, id);
vary_Ii_PINGING_E_MFR5 = tmp_E_MFR_dt(:, id);
vary_Ii_PINGING_I_MFR5 = tmp_I_MFR_dt(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_PINGING_freq = tmp(id, :);
vary_Ie_PINGING_E_MFR = tmp_E_MFR_dt(id, :);
vary_Ie_PINGING_I_MFR = tmp_I_MFR_dt(id, :);

figure(400);hold on;
tmp_E_MFR_dt = reshape(E_MFR_dt, N_y, N_x);
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp_E_MFR_dt);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

figure(401);hold on;
tmp_I_MFR_dt = reshape(I_MFR_dt, N_y, N_x);
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp_I_MFR_dt);
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

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
Perc1 = 'NWCA1_EEEIIEII_IappEx_IappIx_sigmaWNE20_rndV0_tEnd4000';
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0.3_pII0.3_pIIGJ0.004\v4\';
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
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% hSurface = surf(e_freq');
set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);

id = sum(x_lin_correct_unit <= Fixed_Ie);
vary_Ii_PING_freq = tmp(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_PING_freq = tmp(id, :);

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
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0_pII0.3_pIIGJ0.004\v4\';
load(strcat(dir,strcat(Perc, '.mat')))

kappa_dt_for_PING = kappa_dt;
PowerFreq_dt_for_PING = PowerFreq_dt;
Freq_dt_for_PING = Freq_dt;
E_Freq_dt_for_PING = E_Freq_dt;
I_Freq_dt_for_PING = I_Freq_dt;
MFR_dt_for_PING = MFR_dt;

N_x_for_PING = size(unique(x_lin), 1);
N_y_for_PING = size(unique(y_lin), 1);

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
x_lin = unique(x_lin);
y_lin = unique(y_lin);

% x_lin_correct_unit = x_lin;
% y_lin_correct_unit = y_lin;

x_lin_correct_unit = x_lin./nw_A_e.*(1e-12).*(1e+6);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-12).*(1e+6);

N_x = size(x_lin_correct_unit, 1);
N_y = size(y_lin_correct_unit, 1);

b_rhythms = (E_PowerFreq_dt < 0.02) | (isnan(E_PowerFreq_dt));
% b_rhythms = (E_PowerFreq_dt < 10) | (isnan(E_PowerFreq_dt));
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;
I_Freq_dt(b_I_rhythms, 1) = NaN;

b_low_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
b_rhythms = (E_PowerFreq_dt < 10) | (isnan(E_PowerFreq_dt));
E_Freq_dt(b_rhythms & b_low_I_rhythms, 1) = NaN;

figure(4);hold on;
tmp = reshape(I_Freq_dt, N_y, N_x);
% tmp = reshape(E_Freq_dt, N_y, N_x);
tmp_I_MFR_dt = reshape(I_MFR_dt, N_y, N_x);

hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
min_ING_freq = min(min(tmp));
max_ING_freq = max(max(tmp));
min_ING_x = min(min(x_lin_correct_unit));
max_ING_x = max(max(x_lin_correct_unit));
min_ING_y = min(min(y_lin_correct_unit));
max_ING_y = max(max(y_lin_correct_unit));

id = sum(x_lin_correct_unit <= Fixed_Ie);
vary_Ii_ING_freq = tmp(:, id);
vary_Ii_ING_I_MFR = tmp_I_MFR_dt(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_ING_freq = tmp(id, :);

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

figure(4);hold on;
box on
% xlabel('I_{E}');
% ylabel('I_{I}');
axis tight
axis square
view([-37.500000000000000, 30])

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.6 0.8 1 1.2 1.4],'YTickLabel',{'';'';'';'';''});
% set(gca,'ZTick',[30 35 40 45],'ZTickLabel',{'';'';'';''});

% savefig('PINGING_IeIi_Efreq_SigmaI0_5_for_paper2', 'eps');

% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.6 0.8 1 1.2 1.4],'YTickLabel',{'';'';'';'';''});
% set(gca,'ZTick',[10 20 30 40],'ZTickLabel',{'';'';'';''});
% 
% m_savefig('PINGING_IeIi_Efreq_SigmaI0_5_for_paper2_usingEcellsForINGFreq', 'eps');

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

figure(400);hold on;
box on
axis tight
axis square
view([-37.500000000000000, 30])

figure(401);hold on;
box on
axis tight
axis square
view([-37.500000000000000, 30])


%% Create the cross section
% X1 = [Fixed_Ie min_ING_y min_ING_freq];
% X2 = [Fixed_Ie max_ING_y min_ING_freq];
% X3 = [Fixed_Ie max_ING_y max_ING_freq];
% X4 = [Fixed_Ie min_ING_y max_ING_freq];
% makePatches(X1, X2, X3, X4, Patch_FaceAlpha);
% 
% X1 = [min_ING_x Fixed_Ii min_ING_freq];
% X2 = [max_ING_x Fixed_Ii min_ING_freq];
% X3 = [max_ING_x Fixed_Ii max_ING_freq];
% X4 = [min_ING_x Fixed_Ii max_ING_freq];
% makePatches(X1, X2, X3, X4, Patch_FaceAlpha);

figure(100);hold on;
plot(y_lin_correct_unit, vary_Ii_PINGING_freq + 0.1, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
plot(y_lin_correct_unit, vary_Ii_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(y_lin_correct_unit, vary_Ii_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');

% axis tight
axis square
box on
% ylim([min_ING_freq max_ING_freq]);
xlim([0.4981 1.486]);
ylim([30 46]);

set(gca,'XTick',[0.6 0.8 1 1.2 1.4],'XTickLabel',{'';'';'';'';''});
set(gca,'YTick',[30 35 40 45],'YTickLabel',{'';'';'';''});

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% savefig('vary_Ii_freq_SigmaI0_5_for_paper2_Fixed_Ie_2', 'eps');

figure(1000);hold on;
[AX,H1,H2] = plotyy(y_lin_correct_unit, vary_Ii_PINGING_E_MFR, y_lin_correct_unit, vary_Ii_PINGING_I_MFR, 'plot');
box on
% set(get(AX(1),'Ylabel'),'String','MFR_E [Spks/s]') 
% set(get(AX(2),'Ylabel'),'String','MFR_I [Spks/s]') 
set(get(AX(1),'Ylabel'),'String','') 
set(get(AX(2),'Ylabel'),'String','') 

set(AX(1), 'xlim', [0.4981 1.486]);
set(AX(2), 'xlim', [0.4981 1.486]);
axis(AX(1), 'square');
axis(AX(2), 'square');

set(H1, 'LineStyle', 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.609 0], 'MarkerEdgeColor', [0 0.609 0]);
set(H2, 'LineStyle', 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.6 1.0 0.6], 'MarkerEdgeColor', [0.6 1.0 0.6]);

set(AX(1),'XTick',[0.6 0.8 1 1.2 1.4],'XTickLabel',{'';'';'';'';''});
set(AX(1),'YTick',[0 1 2 3 4 5],'YTickLabel',{'';'';'';'';'';''});

set(AX(2),'XTick',[0.6 0.8 1 1.2 1.4],'XTickLabel',{'';'';'';'';''});
set(AX(2),'YTick',[34 36 38 40 42 44],'YTickLabel',{'';'';'';'';'';''});

make_me_pretty(gcf, ...
    AX(1), 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

make_me_pretty(gcf, ...
    AX(2), 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% m_savefig('vary_Ii_freq_SigmaI0_5_for_paper2_Fixed_Ie_2_MFR_E_MFR_I', 'eps');

figure(101);hold on;
plot(x_lin_correct_unit, vary_Ie_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
plot(x_lin_correct_unit, vary_Ie_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(x_lin_correct_unit, vary_Ie_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');

axis tight
axis square
box on
% ylim([min_ING_freq max_ING_freq]);
ylim([30 46]);

set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
set(gca,'YTick',[30 35 40 45],'YTickLabel',{'';'';'';''});

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

Fixed_Ie = 2.0;
Fixed_Ii = 0.85;
% savefig('vary_Ie_freq_SigmaI0_5_for_paper2_Fixed_Ii_0_85', 'eps');

figure(1010);hold on;
[AX,H1,H2] = plotyy(x_lin_correct_unit, vary_Ie_PINGING_E_MFR, x_lin_correct_unit, vary_Ie_PINGING_I_MFR, 'plot');
box on
% set(get(AX(1),'Ylabel'),'String','MFR_E [Spks/s]') 
% set(get(AX(2),'Ylabel'),'String','MFR_I [Spks/s]') 
set(get(AX(1),'Ylabel'),'String','') 
set(get(AX(2),'Ylabel'),'String','') 

set(AX(1), 'xlim', [0 2.316]);
set(AX(2), 'xlim', [0 2.316]);
axis(AX(1), 'square');
axis(AX(2), 'square');

set(H1, 'LineStyle', 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.609 0], 'MarkerEdgeColor', [0 0.609 0]);
set(H2, 'LineStyle', 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.6 1.0 0.6], 'MarkerEdgeColor', [0.6 1.0 0.6]);

set(AX(1),'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
set(AX(1),'YTick',[0 1 2 3 4 5],'YTickLabel',{'';'';'';'';'';''});

set(AX(2),'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
set(AX(2),'YTick',[20 24 28 32 36 40],'YTickLabel',{'';'';'';'';'';''});


make_me_pretty(gcf, ...
    AX(1), 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

make_me_pretty(gcf, ...
    AX(2), 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% m_savefig('vary_Ie_freq_SigmaI0_5_for_paper2_Fixed_Ii_0_85_MFR_E_MFR_I', 'eps');

y_lin_correct_unit = y_lin_correct_unit.*wb_A_i./(1e-12)./(1e+6);

figure(1012);hold on;
vary_Ii_ING_freq = tmp(:, id);
vary_Ii_ING_I_MFR = tmp_I_MFR_dt(:, id);
plot(y_lin_correct_unit, vary_Ii_ING_I_MFR./(vary_Ii_ING_freq), 'r*');


figure(1011);hold on;
% plot(y_lin_correct_unit, , 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR./(vary_Ii_PINGING_freq), 'r*');

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

axis square
grid on
box on

% xlim([0.4981 1.5]);

maximize_a_fig(gcf);
% m_savefig('type1INT_ratioFreqMFRI_for_1reviewer', 'eps');

figure(2012);
y_lin_correct_unit = y_lin_correct_unit./wb_A_i.*(1e-12).*(1e+6);
subplot(1, 2, 1);hold on;
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR1, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR2, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR3, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR4, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_I_MFR5, 'b');
grid on

subplot(1, 2, 2);hold on;
plot(y_lin_correct_unit, vary_Ii_PINGING_E_MFR1, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_E_MFR2, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_E_MFR3, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_E_MFR4, 'r');
plot(y_lin_correct_unit, vary_Ii_PINGING_E_MFR5, 'b');
grid on

end

