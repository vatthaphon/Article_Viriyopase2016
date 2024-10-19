function my_plot
clc;clear all;
close all;
format long

Perc = 'X';
PINGING_FaceAlpha = 0.5 + 0.0;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.0;
Ii_min = 90;
Ii_max = 270;

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

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;

figure(4);hold on;
% tmp = reshape(Freq_dt, N_y, N_x);
tmp = reshape(E_Freq_dt, N_y, N_x);
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

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
running_i = 1;
for i = 1:N_x_for_PING
    for j = 1:N_y_for_PING
        E_kappa_dtMAT(i,j) = kappa_dt_for_PING(running_i, 1);
        I_kappa_dtMAT(i,j) = kappa_dt_for_PING(running_i, 2);
        PowerFreq_dtMAT(i,j) = PowerFreq_dt_for_PING(running_i, 1);
        Freq_dtMAT(i,j) = Freq_dt_for_PING(running_i, 1);
        E_Freq_dtMAT(i,j) = E_Freq_dt_for_PING(running_i, 1);
        I_Freq_dtMAT(i,j) = I_Freq_dt_for_PING(running_i, 1);
        E_MFR_dtMAT(i,j) = MFR_dt_for_PING(running_i, 1);
        I_MFR_dtMAT(i,j) = MFR_dt_for_PING(running_i, 2);
        
        running_i = running_i + 1;
    end
end

for i = 1:1:N_x
    for j = 1:1:N_y
        freq(i, j) = Freq_dtMAT(i, 1);
        e_freq(i, j) = E_Freq_dtMAT(i, 1);
        e_mfr(i, j) = E_MFR_dtMAT(i, 1);
        i_mfr(i, j) = I_MFR_dtMAT(i, 1);
        e_kappa(i, j) = E_kappa_dtMAT(i, 1);
        i_kappa(i, j) = I_kappa_dtMAT(i, 1);        
    end
end

figure(4);hold on;
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, e_freq');
% hSurface = surf(e_freq');
set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);

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

b_rhythms = (E_PowerFreq_dt < 0.003) | (isnan(E_PowerFreq_dt));
E_PowerFreq_dt(b_rhythms, 1) = NaN;
E_Freq_dt(b_rhythms, 1) = NaN;
MFR_dt(b_rhythms, 1) = NaN;
kappa_dt(b_rhythms, 1) = NaN;

b_I_rhythms = (I_PowerFreq_dt < 1) | (isnan(I_PowerFreq_dt));
MFR_dt(b_I_rhythms, 2) = NaN;
kappa_dt(b_I_rhythms, 2) = NaN;

figure(4);hold on;
tmp = reshape(I_Freq_dt, N_y, N_x);
hSurface = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

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
% set(gca,'YTick',[0.6 0.8 1 1.2],'YTickLabel',{'';'';'';''});
% set(gca,'ZTick',[25 30 35 40],'ZTickLabel',{'';'';'';''});

% savefig('PINGING_freq', 'eps');

%savefig('PINGING_IeIi_Efreq', 'eps');

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
end


