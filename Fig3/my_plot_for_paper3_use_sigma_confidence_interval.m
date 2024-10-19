function my_plot_for_paper3_use_sigma_confidence_interval
clc;clear all;close all;
format long

N_rndV = 1;

for rndV = 0:1:(N_rndV - 1)
    
    Perc = strcat('X');
    dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\bw_hh\EEIEII\NWCA1_gIIxi_gIIGJx_0.01_20_41_IappI1280_IappE1258.9\v0\data\';
    Asyn_th = 40;
    
    load(strcat(dir,strcat(Perc, '.mat')))
    
    tmp_I_PowerFreq_dt = reshape(I_PowerFreq_dt(:, 1), N_y, N_x)';
    tmp_I_Freq_dt = reshape(I_Freq_dt(:, 1), N_y, N_x)';
    
    Rhythm_certainty_Stan_Raoul_power = tmp_I_PowerFreq_dt;
    tmp_Rhythm_certainty_Stan_Raoul_power = Rhythm_certainty_Stan_Raoul_power ;
    tmp_Freq_dt = tmp_I_Freq_dt;
    
end

power_for_asyn_state = Rhythm_certainty_Stan_Raoul_power(Rhythm_certainty_Stan_Raoul_power < Asyn_th);

tmp_mean = mean(power_for_asyn_state);
tmp_sigma = std(power_for_asyn_state);

tmp_Freq_dt_1sigma = tmp_Freq_dt;
tmp_Freq_dt_2sigma = tmp_Freq_dt;
tmp_Freq_dt_3sigma = tmp_Freq_dt;
tmp_Freq_dt_4sigma = tmp_Freq_dt;
tmp_Freq_dt_5sigma = tmp_Freq_dt;
tmp_Freq_dt_6sigma = tmp_Freq_dt;
tmp_Freq_dt_7sigma = tmp_Freq_dt;
tmp_Freq_dt_8sigma = tmp_Freq_dt;

tmp_Freq_dt_1sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 1*tmp_sigma)) = NaN;
tmp_Freq_dt_2sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 2*tmp_sigma)) = NaN;
tmp_Freq_dt_3sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 3*tmp_sigma)) = NaN;
tmp_Freq_dt_4sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 4*tmp_sigma)) = NaN;
tmp_Freq_dt_5sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 5*tmp_sigma)) = NaN;
tmp_Freq_dt_6sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 6*tmp_sigma)) = NaN;
tmp_Freq_dt_7sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 7*tmp_sigma)) = NaN;
tmp_Freq_dt_8sigma(tmp_Rhythm_certainty_Stan_Raoul_power < (tmp_mean + 8*tmp_sigma)) = NaN;

wb_A_i = (18069*1e-8); % [cm^2]
x_lin = 0.01*((20.0/0.01).^(x_lin/(41 - 1.0)));

x_lin = unique(x_lin);
y_lin = unique(y_lin);

x_lin_correct_unit = x_lin./wb_A_i.*(1e-9).*(1e+3);
y_lin_correct_unit = y_lin./wb_A_i.*(1e-9).*(1e+3);

% x_lin_correct_unit = x_lin;
% y_lin_correct_unit = y_lin;

my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, tmp_Freq_dt_3sigma, ...
    'is_XLog', 1, 'is_YLog', 0, ...
    'CLimBegin', 50, 'CLimEnd', 70, ...
    'is_ShowEdge', 0);
colorbar

make_me_pretty(gcf, ...
    gca, 15*1.5, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
axis square

m_savefig('gII_vs_gGJII_3sigma', 'eps');

end


