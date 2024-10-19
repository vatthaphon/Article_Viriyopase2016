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
Fixed_Ie = 1.25;
Fixed_Ii = 0.95;

nw_A_e = (21590*1e-8); % [cm^2]
wb_A_i = (18069*1e-8); % [cm^2]

%% PING + ING
dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0.3_pII0.3_pIIGJ0.004\v5\';
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
t_lags_dt = t_lags_dt(good_y_lin, 1);

t_lags_dt(t_lags_dt>=10) = NaN;
t_lags_dt(t_lags_dt<=-10) = NaN;

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

% Create a 3D matrix
x_lin = unique(x_lin);
y_lin = unique(y_lin);

x_lin_correct_unit = x_lin;
y_lin_correct_unit = y_lin;

x_lin_correct_unit = x_lin./nw_A_e.*(1e-12).*(1e+6);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-12).*(1e+6);

N_x = size(x_lin_correct_unit, 1);
N_y = size(y_lin_correct_unit, 1);

b_rhythms = (E_PowerFreq_dt < 0.01) | (isnan(E_PowerFreq_dt)) | (isnan(MFR_dt(:, 1)))  | (MFR_dt(:, 1) <= 0.15) ;
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;
t_lags_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;

figure(4);hold on;
% tmp = reshape(Freq_dt, N_y, N_x);
tmp = reshape(E_Freq_dt, N_y, N_x);
% tmp = reshape(t_lags_dt, N_y, N_x);
t_lags_dt_PINGING = tmp;
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

id = sum(x_lin_correct_unit <= Fixed_Ie);
vary_Ii_PINGING_freq = tmp(:, id);

id = sum(y_lin_correct_unit <= Fixed_Ii);
vary_Ie_PINGING_freq = tmp(id, :);

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
dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0.3_pII0.3_pIIGJ0.004\v4\';
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
    t_lags_dtMAT(i,1) = t_lags_dt(running_i, 1);
    
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
        t_lags(i, j) = t_lags_dtMAT(i, 1);
    end
end

figure(4);hold on;
tmp = e_freq';
% tmp = t_lags';
t_lags_dt_PING = tmp;
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
dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_IappIx_sigmaWNE20_pEI0_pII0.3_pIIGJ0.004\v3\';
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
t_lags_dt = t_lags_dt(good_y_lin, 1);

t_lags_dt(t_lags_dt >= 10) = NaN;
t_lags_dt(t_lags_dt <= -10) = NaN;

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

% Create a 3D matrix
x_lin = unique(x_lin);
y_lin = unique(y_lin);

x_lin_correct_unit = x_lin;
y_lin_correct_unit = y_lin;

x_lin_correct_unit = x_lin./nw_A_e.*(1e-12).*(1e+6);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-12).*(1e+6);

N_x = size(x_lin_correct_unit, 1);
N_y = size(y_lin_correct_unit, 1);

b_rhythms = (E_PowerFreq_dt < 0.003) | (isnan(E_PowerFreq_dt) | (MFR_dt(:, 1) > E_Freq_dt(:, 1)));
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;
t_lags_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;

figure(4);hold on;
tmp = reshape(I_Freq_dt, N_y, N_x);
% tmp = reshape(t_lags_dt, N_y, N_x);
t_lags_dt_ING = tmp;
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
% set(gca,'ZTick',[25 30 35 40],'ZTickLabel',{'';'';'';''});

% savefig('PINGING_freq', 'eps');

% savefig('PINGING_IeIi_Efreq', 'eps');

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

% figure(100);hold on;
% plot(y_lin_correct_unit, vary_Ii_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
% plot(y_lin_correct_unit, vary_Ii_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
% plot(y_lin_correct_unit, vary_Ii_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
% 
% axis tight
% axis square
% box on
% ylim([min_ING_freq max_ING_freq]);
% 
% set(gca,'XTick',[0.6 0.8 1 1.2 1.4],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[25 30 35 40],'YTickLabel',{'';'';'';''});
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
% % savefig('vary_Ii_freq', 'eps');
% 
% figure(101);hold on;
% plot(x_lin_correct_unit, vary_Ie_PINGING_freq, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g');
% plot(x_lin_correct_unit, vary_Ie_PING_freq, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
% plot(x_lin_correct_unit, vary_Ie_ING_freq, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
% 
% axis tight
% axis square
% box on
% ylim([min_ING_freq max_ING_freq]);
% 
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[25 30 35 40],'YTickLabel',{'';'';'';''});
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
% % savefig('vary_Ie_freq', 'eps');

%% Calculate the tendency.
t_lags_dt_PING;
t_lags_dt_PINGING;
t_lags_dt_ING;

[N_y, N_x] = size(t_lags_dt_PING);

% Avarage along the x-axis
t_lags_dt_PING_avg_along_x = NaN(1, N_x);
for i = 1:1:N_x
    n = 0;
    for j = 1:1:N_y
        if (isnan(t_lags_dt_PING(j, i)) == 0)
            if (isnan(t_lags_dt_PING_avg_along_x(1, i)) == 1)
                t_lags_dt_PING_avg_along_x(1, i) = t_lags_dt_PING(j, i);            
            else
                t_lags_dt_PING_avg_along_x(1, i) = t_lags_dt_PING_avg_along_x(1, i) + t_lags_dt_PING(j, i);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_PING_avg_along_x(1, i) = t_lags_dt_PING_avg_along_x(1, i)/n;            
end

figure(99);hold on;
plot(x_lin_correct_unit, t_lags_dt_PING_avg_along_x, 'r*')

% Avarage along the y-axis
t_lags_dt_PING_avg_along_y = NaN(1, N_y);
for i = 1:1:N_y
    n = 0;
    for j = 1:1:N_x
        if (isnan(t_lags_dt_PING(i, j)) == 0)
            if (isnan(t_lags_dt_PING_avg_along_y(1, i)) == 1)
                t_lags_dt_PING_avg_along_y(1, i) = t_lags_dt_PING(i, j);            
            else
                t_lags_dt_PING_avg_along_y(1, i) = t_lags_dt_PING_avg_along_y(1, i) + t_lags_dt_PING(i, j);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_PING_avg_along_y(1, i) = t_lags_dt_PING_avg_along_y(1, i)/n;            
end

figure(999);hold on;
plot(y_lin_correct_unit, t_lags_dt_PING_avg_along_y, 'r*')

[N_y, N_x] = size(t_lags_dt_ING);

% Avarage along the x-axis
t_lags_dt_ING_avg_along_x = NaN(1, N_x);
for i = 1:1:N_x
    n = 0;
    for j = 1:1:N_y
        if (isnan(t_lags_dt_ING(j, i)) == 0)
            if (isnan(t_lags_dt_ING_avg_along_x(1, i)) == 1)
                t_lags_dt_ING_avg_along_x(1, i) = t_lags_dt_ING(j, i);            
            else
                t_lags_dt_ING_avg_along_x(1, i) = t_lags_dt_ING_avg_along_x(1, i) + t_lags_dt_ING(j, i);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_ING_avg_along_x(1, i) = t_lags_dt_ING_avg_along_x(1, i)/n;            
end

figure(99);hold on;
plot(x_lin_correct_unit, t_lags_dt_ING_avg_along_x, 'b*')

% Avarage along the y-axis
t_lags_dt_ING_avg_along_y = NaN(1, N_y);
for i = 1:1:N_y
    n = 0;
    for j = 1:1:N_x
        if (isnan(t_lags_dt_ING(i, j)) == 0)
            if (isnan(t_lags_dt_ING_avg_along_y(1, i)) == 1)
                t_lags_dt_ING_avg_along_y(1, i) = t_lags_dt_ING(i, j);            
            else
                t_lags_dt_ING_avg_along_y(1, i) = t_lags_dt_ING_avg_along_y(1, i) + t_lags_dt_ING(i, j);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_ING_avg_along_y(1, i) = t_lags_dt_ING_avg_along_y(1, i)/n;            
end

figure(999);hold on;
plot(y_lin_correct_unit, t_lags_dt_ING_avg_along_y, 'b*')

[N_y, N_x] = size(t_lags_dt_PINGING);

% Avarage along the x-axis
t_lags_dt_PINGING_avg_along_x = NaN(1, N_x);
for i = 1:1:N_x
    n = 0;
    for j = 1:1:N_y
        if (isnan(t_lags_dt_PINGING(j, i)) == 0)
            if (isnan(t_lags_dt_PINGING_avg_along_x(1, i)) == 1)
                t_lags_dt_PINGING_avg_along_x(1, i) = t_lags_dt_PINGING(j, i);            
            else
                t_lags_dt_PINGING_avg_along_x(1, i) = t_lags_dt_PINGING_avg_along_x(1, i) + t_lags_dt_PINGING(j, i);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_PINGING_avg_along_x(1, i) = t_lags_dt_PINGING_avg_along_x(1, i)/n;            
end

figure(99);hold on;
plot(x_lin_correct_unit, t_lags_dt_PINGING_avg_along_x, 'g*')

% Avarage along the y-axis
t_lags_dt_PINGING_avg_along_y = NaN(1, N_y);
for i = 1:1:N_y
    n = 0;
    for j = 1:1:N_x
        if (isnan(t_lags_dt_PINGING(i, j)) == 0)
            if (isnan(t_lags_dt_PINGING_avg_along_y(1, i)) == 1)
                t_lags_dt_PINGING_avg_along_y(1, i) = t_lags_dt_PINGING(i, j);            
            else
                t_lags_dt_PINGING_avg_along_y(1, i) = t_lags_dt_PINGING_avg_along_y(1, i) + t_lags_dt_PINGING(i, j);            
            end
            n = n + 1;
        end
    end
    t_lags_dt_PINGING_avg_along_y(1, i) = t_lags_dt_PINGING_avg_along_y(1, i)/n;            
end

figure(999);hold on;
plot(y_lin_correct_unit, t_lags_dt_PINGING_avg_along_y, 'g*')

figure(99)
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

box on
ylabel('E-cells lead I-cells by [ms]')

% axis tight
% savefig('vary_Ii_freq', 'eps');

figure(999)
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

box on
xlabel('Current to interneurons [microA/cm^2]');
ylabel('Frequency [Hz]');
% ylabel('E-cells lead I-cells by [ms]');
axis tight
legend('PING','ING','PING+ING');
% savefig('vary_Ii_freq', 'eps');

end

function makePatches(X1, X2, X3, X4, Patch_FaceAlpha)
%   INPUT:
%     X1   - 1 dimensional 1-by-3, (x, y, z)
%     X2   - 1 dimensional 1-by-3, (x, y, z)
%     X3   - 1 dimensional 1-by-3, (x, y, z)
%     X4   - 1 dimensional 1-by-3, (x, y, z)

patch('vertices', [X1(1, 1) X1(1, 2) X1(1, 3); X2(1, 1) X2(1, 2) X2(1, 3); X3(1, 1) X3(1, 2) X3(1, 3); X4(1, 1) X4(1, 2) X4(1, 3)], ...
    'faces',[1 2 3 4], 'facecolor',[.5 .5 .5],'FaceAlpha', Patch_FaceAlpha);

end
