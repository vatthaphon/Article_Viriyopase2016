function combine
clc;clear all;close all;
format long

tmp_N_x = 41;
tmp_N_y = 41;

tmp_E_Freq_dt = [];
tmp_I_Freq_dt = [];
tmp_Freq_dt = [];
tmp_MFR_dt_E = [];
tmp_MFR_dt_I = [];
tmp_E_PowerFreq_dt = [];
tmp_I_PowerFreq_dt = [];
tmp_PowerFreq_dt = [];
tmp_t_lags_dt = [];
tmp_E_Rhythm_certainty = [];
tmp_I_Rhythm_certainty = [];
tmp_Rhythm_certainty = [];
tmp_kappa_dt_E = [];
tmp_kappa_dt_I = [];
tmp_x_lin = linspace(0.0, 500.0, tmp_N_x);
tmp_y_lin = linspace(0, 20, tmp_N_y);

i = 1;
for X = tmp_x_lin
    
Perc = strcat('X_', num2str(X));

dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_sigmaWNEx_pII0.3_pIIGJ0.004\v2\';
load(strcat(dir,strcat(Perc, '.mat')))

tmp_Freq_dt             = [tmp_Freq_dt; Freq_dt(1:tmp_N_y, 1)];
tmp_MFR_dt_E            = [tmp_MFR_dt_E; MFR_dt(1:tmp_N_y, 1)];
tmp_MFR_dt_I            = [tmp_MFR_dt_I; MFR_dt(1:tmp_N_y, 2)];
tmp_PowerFreq_dt        = [tmp_PowerFreq_dt; PowerFreq_dt(1:tmp_N_y, 1)];
tmp_t_lags_dt           = [tmp_t_lags_dt; t_lags_dt(1:tmp_N_y, 1)];
tmp_Rhythm_certainty    = [tmp_Rhythm_certainty; Rhythm_certainty(1:tmp_N_y, 1)];
tmp_kappa_dt_E          = [tmp_kappa_dt_E; kappa_dt(1:tmp_N_y, 1)];
tmp_kappa_dt_I          = [tmp_kappa_dt_I; kappa_dt(1:tmp_N_y, 2)];
tmp_E_Freq_dt           = [tmp_E_Freq_dt; E_Freq_dt(1:tmp_N_y, 1)];
tmp_I_Freq_dt           = [tmp_I_Freq_dt; I_Freq_dt(1:tmp_N_y, 1)];
tmp_E_PowerFreq_dt      = [tmp_E_PowerFreq_dt; E_PowerFreq_dt(1:tmp_N_y, 1)];
tmp_I_PowerFreq_dt      = [tmp_I_PowerFreq_dt; I_PowerFreq_dt(1:tmp_N_y, 1)];
tmp_E_Rhythm_certainty  = [tmp_E_Rhythm_certainty; E_Rhythm_certainty(1:tmp_N_y, 1)];
tmp_I_Rhythm_certainty  = [tmp_I_Rhythm_certainty; I_Rhythm_certainty(1:tmp_N_y, 1)];

i = i + 1;
end

tmp_x_lin = [];
tmp_y_lin = [];

for X = linspace(0.0, 500.0, tmp_N_x)
for Y = linspace(0, 20, tmp_N_y)
    tmp_x_lin = [tmp_x_lin; X];
    tmp_y_lin = [tmp_y_lin; Y];
end
end

N_x = tmp_N_x;
N_y = tmp_N_y;
x_lin = tmp_x_lin;
y_lin = tmp_y_lin;
Freq_dt = tmp_Freq_dt;
t_lags_dt = tmp_t_lags_dt;
MFR_dt = [tmp_MFR_dt_E tmp_MFR_dt_I];
PowerFreq_dt = tmp_PowerFreq_dt;
Rhythm_certainty = tmp_Rhythm_certainty;
kappa_dt = [tmp_kappa_dt_E tmp_kappa_dt_I];
E_Freq_dt = tmp_E_Freq_dt;
I_Freq_dt = tmp_I_Freq_dt;
E_PowerFreq_dt = tmp_E_PowerFreq_dt;
I_PowerFreq_dt = tmp_I_PowerFreq_dt;
E_Rhythm_certainty = tmp_E_Rhythm_certainty;
I_Rhythm_certainty  = tmp_I_Rhythm_certainty;

save('X.mat', 'x_lin', 'y_lin', 'E_Freq_dt', 'I_Freq_dt', 'E_PowerFreq_dt', 'I_PowerFreq_dt', 'E_Rhythm_certainty', 'I_Rhythm_certainty', 't_lags_dt', 'Rhythm_certainty', 'MFR_dt', 'kappa_dt', 'PowerFreq_dt', 'Freq_dt', 'N_x', 'N_y');

end


