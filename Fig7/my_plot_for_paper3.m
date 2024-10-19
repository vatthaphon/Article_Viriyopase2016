function my_plot
clc;clear all;
close all;
format long

nw_A_e = (21590*1e-8); % [cm^2]
wb_A_i = (18069*1e-8); % [cm^2]
factor1 = (1e-12);
factor2 = (1e+6);

Perc = 'X';
Patch_FaceAlpha = 0.4;
PINGING_FaceAlpha = 0.5 + 0.2;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.4;

Ii_min = 0.0*wb_A_i/((1e-12).*(1e+6));
Ii_max = 20.0*wb_A_i/((1e-12).*(1e+6));

Ie_min = 6.0;
% Ie_min = 0.0;
Ie_max = 30.0;

Fixed_Ii = 7.5*wb_A_i/((1e-12).*(1e+6));
Fixed_Ie = 17*nw_A_e/((1e-12).*(1e+6));

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

good_y_lin = (Ii_min <= y_lin) & (y_lin <= Ii_max) ;
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

% Create a 3D matrix
x_lin = begin_Val13*((end_Val13/begin_Val13).^(x_lin/(len_Val13 - 1.0)));

x_lin_correct_unit = x_lin./nw_A_e.*factor1.*factor2;
y_lin_correct_unit = y_lin./wb_A_i.*factor1.*factor2;

good_x_lin = (Ie_min <= x_lin_correct_unit) & (x_lin_correct_unit <= Ie_max) ;
x_lin_correct_unit = x_lin_correct_unit(good_x_lin);
y_lin_correct_unit = y_lin_correct_unit(good_x_lin);
E_Freq_dt = E_Freq_dt(good_x_lin, 1);
I_Freq_dt = I_Freq_dt(good_x_lin, 1);
E_PowerFreq_dt = E_PowerFreq_dt(good_x_lin, 1);
I_PowerFreq_dt = I_PowerFreq_dt(good_x_lin, 1);
E_MFR_dt = MFR_dt(good_x_lin, 1);
I_MFR_dt = MFR_dt(good_x_lin, 2);
E_kappa_dt = kappa_dt(good_x_lin, 1);
I_kappa_dt = kappa_dt(good_x_lin, 2);

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

x_lin_correct_unit = unique(x_lin_correct_unit);  % Ie
y_lin_correct_unit = unique(y_lin_correct_unit);  % Ii

N_x_for_PING = size(unique(x_lin_correct_unit), 1);
N_y_for_PING = size(unique(y_lin_correct_unit), 1);

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
tmp_MFR_E = reshape(MFR_dt(:, 1), N_y, N_x);
tmp_MFR_I = reshape(MFR_dt(:, 2), N_y, N_x);

% z_const2 = tmp(41, :);
% max_y = max(y_lin_correct_unit);
% plot3(x_lin_correct_unit, ones(size(x_lin_correct_unit, 1), 1)*max_y, z_const2, 'LineWidth', 20, 'Color', [0 0 0]);
hSurface1 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

% id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
% vary_Ii_PINGING_freq = tmp(:, id);
% vary_Ie_PINGING_MFR_E = tmp_MFR_E(:, id);
% vary_Ie_PINGING_MFR_I = tmp_MFR_I(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_PINGING_freq = tmp(id, :);
% vary_Ie_PINGING_MFR_E = tmp_MFR_E(id, :);
% vary_Ie_PINGING_MFR_I = tmp_MFR_I(id, :);

% figure(5);hold on;
% tmp = reshape(E_MFR_dt, N_y, N_x);
% Fig5_PINGING_x_lin_correct_unit = x_lin_correct_unit;
% Fig5_PINGING_y_lin_correct_unit = y_lin_correct_unit;
% Fig5_PINGING_tmp = tmp;

% figure(6);hold on;
% tmp = reshape(I_MFR_dt, N_y, N_x);
% hSurface3 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

%% PING
Perc1 = 'X';
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEEIIEII\NWCA1_IappExi_10_10000_41_IappI0_sigmaWNE60_sigmaWNI0.5\v0\data\';
load(strcat(dir, strcat(Perc1, '.mat')))

running_i = 1;
for i = 1:N_x
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

for i = 1:1:N_x
    for j = 1:1:size(y_lin_correct_unit, 1)
        freq(i, j) = Freq_dtMAT(i, 1);
        e_freq(i, j) = E_Freq_dtMAT(i, 1);
        e_mfr(i, j) = E_MFR_dtMAT(i, 1);
        i_mfr(i, j) = I_MFR_dtMAT(i, 1);
        e_kappa(i, j) = E_kappa_dtMAT(i, 1);
        i_kappa(i, j) = I_kappa_dtMAT(i, 1);        
    end
end

x_lin = begin_Val13*((end_Val13/begin_Val13).^(x_lin/(len_Val13 - 1.0)));
x_lin_correct_unit = x_lin./nw_A_e.*factor1.*factor2;

good_x_lin = (Ie_min <= x_lin_correct_unit) & (x_lin_correct_unit <= Ie_max) ;
x_lin_correct_unit = x_lin_correct_unit(good_x_lin);

x_lin_correct_unit = log10(x_lin_correct_unit);

figure(4);hold on;
tmp = e_freq';

tmp = tmp(:, good_x_lin);
hSurface4 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

% id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
% vary_Ii_PING_freq = tmp(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_PING_freq = tmp(id, :);
% 
% figure(5);hold on;
% tmp = e_mfr';
% Fig5_PING_x_lin_correct_unit = x_lin_correct_unit;
% Fig5_PING_y_lin_correct_unit = y_lin_correct_unit;
% Fig5_PING_tmp = tmp;
% 
% figure(6);hold on;
% tmp = i_mfr';
% hSurface6 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

%% ING
dir = 'C:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEIEII\NWCA1_IappExi_IappIx_10_10000_41_sigmaWNE60_sigmaWNI0.5\v0\data\';
load(strcat(dir,strcat(Perc, '.mat')))

good_y_lin = (Ii_min <= y_lin) & (y_lin <= Ii_max) ;
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

x_lin = begin_Val13*((end_Val13/begin_Val13).^(x_lin/(len_Val13 - 1.0)));

x_lin_correct_unit = x_lin./nw_A_e.*factor1.*factor2;
y_lin_correct_unit = y_lin./wb_A_i.*factor1.*factor2;

good_x_lin = (Ie_min <= x_lin_correct_unit) & (x_lin_correct_unit <= Ie_max) ;
x_lin_correct_unit = x_lin_correct_unit(good_x_lin);
y_lin_correct_unit = y_lin_correct_unit(good_x_lin);
E_Freq_dt = E_Freq_dt(good_x_lin, 1);
I_Freq_dt = I_Freq_dt(good_x_lin, 1);
E_PowerFreq_dt = E_PowerFreq_dt(good_x_lin, 1);
I_PowerFreq_dt = I_PowerFreq_dt(good_x_lin, 1);
E_MFR_dt = MFR_dt(good_x_lin, 1);
I_MFR_dt = MFR_dt(good_x_lin, 2);
E_kappa_dt = kappa_dt(good_x_lin, 1);
I_kappa_dt = kappa_dt(good_x_lin, 2);

MFR_dt = [E_MFR_dt I_MFR_dt];
kappa_dt = [E_kappa_dt I_kappa_dt];

x_lin_correct_unit = unique(x_lin_correct_unit);  % Ie
y_lin_correct_unit = unique(y_lin_correct_unit);  % Ii

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

x_lin_correct_unit = log10(x_lin_correct_unit);

figure(4);hold on;
% tmp = reshape(I_Freq_dt, N_y, N_x);
tmp = reshape(E_Freq_dt, N_y, N_x);
hSurface7 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);
% 
% id = sum(x_lin_correct_unit_tmp <= Fixed_Ie);
% vary_Ii_ING_freq = tmp(:, id);
% 
% id = sum(y_lin_correct_unit <= Fixed_Ii);
% vary_Ie_ING_freq = tmp(id, :);
% 
% figure(5);hold on;
% tmp = reshape(MFR_dt(:, 1), N_y, N_x);
% Fig5_ING_x_lin_correct_unit = x_lin_correct_unit;
% Fig5_ING_y_lin_correct_unit = y_lin_correct_unit;
% Fig5_ING_tmp = tmp;
% 
% figure(6);hold on;
% tmp = reshape(MFR_dt(:, 2), N_y, N_x);
% hSurface9 = surf(x_lin_correct_unit, y_lin_correct_unit, tmp);

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

set(hSurface1,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);
set(hSurface4,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);
set(hSurface7,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

figure(4)
set(gca, 'YTick', [6.6 6.8 7 7.2 7.4 7.6], 'YTickLabel',{'';'';'';'';'';''});
% set(gca, 'ZTick', [56 60 64 68], 'ZTickLabel',{'';'';'';''});
maximize_a_fig(gcf);
% m_savefig('PINGING_IeIi_Efreq_for_paper3', 'eps');

end

