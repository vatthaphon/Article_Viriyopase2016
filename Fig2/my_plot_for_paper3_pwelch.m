function my_plot_for_paper3_pwelch
clc;clear all;close all;
format long

% N_rndV = 10;
N_rndV = 1;

first_avg_pwd = zeros(26, 41);
first_avg_freq = zeros(26, 41);
first_avg_MFR_I = zeros(26, 41);

second_avg_pwd = zeros(15, 41);
second_avg_freq = zeros(15, 41);
second_avg_MFR_I = zeros(15, 41);

for rndV = 1
    Perc = strcat('X_rndV', num2str(rndV), '_pwelch');
    
    dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data\';
    load(strcat(dir,strcat(Perc, '.mat')))
    
%     tmp_PowerFreq_dt = reshape(PowerFreq_dt(:, 1), N_y, N_x)';
%     tmp_Freq_dt = reshape(Freq_dt(:, 1), N_y, N_x)';
    tmp_PowerFreq_dt = reshape(I_PowerFreq_dt(:, 1), N_y, N_x)';
    tmp_Freq_dt = reshape(I_Freq_dt(:, 1), N_y, N_x)';
    
    tmp_MFR_dt_I_cells = reshape(MFR_dt(:, 2), N_y, N_x)';
    
    first_tmp_PowerFreq_dt = tmp_PowerFreq_dt(1:26, :);
    first_tmp_Freq_dt = tmp_Freq_dt(1:26, :);
    first_tmp_MFR_dt_I_cells = tmp_MFR_dt_I_cells(1:26, :);
    
    second_tmp_PowerFreq_dt = tmp_PowerFreq_dt(27:41, :);    
    second_tmp_Freq_dt = tmp_Freq_dt(27:41, :);    
    second_tmp_MFR_dt_I_cells = tmp_MFR_dt_I_cells(27:41, :);

%     bad_id = (20 > first_tmp_PowerFreq_dt);    
bad_id = (0 > first_tmp_PowerFreq_dt);    
    first_tmp_PowerFreq_dt(bad_id) = NaN;
    first_tmp_Freq_dt(bad_id) = NaN;
    first_tmp_MFR_dt_I_cells(bad_id) = NaN;
    
%     bad_id = (25 > second_tmp_PowerFreq_dt);    
bad_id = (0 > second_tmp_PowerFreq_dt);    
    second_tmp_PowerFreq_dt(bad_id) = NaN;
    second_tmp_Freq_dt(bad_id) = NaN;
    second_tmp_MFR_dt_I_cells(bad_id) = NaN;
    
    first_avg_pwd = first_avg_pwd + first_tmp_PowerFreq_dt;
    first_avg_freq = first_avg_freq + first_tmp_Freq_dt;
    first_avg_MFR_I = first_avg_MFR_I + first_tmp_MFR_dt_I_cells;
    
    second_avg_pwd = second_avg_pwd + second_tmp_PowerFreq_dt;
    second_avg_freq = second_avg_freq + second_tmp_Freq_dt;
    second_avg_MFR_I = second_avg_MFR_I + second_tmp_MFR_dt_I_cells;    
end

avg_pwd = [first_avg_pwd;second_avg_pwd];
avg_freq = [first_avg_freq;second_avg_freq];
avg_MFR_I = [first_avg_MFR_I;second_avg_MFR_I];

avg_pwd = avg_pwd/N_rndV;
avg_freq = avg_freq/N_rndV;
avg_MFR_I = avg_MFR_I/N_rndV;

% wb_A_i = (18069*1e-8); % [cm^2]
% x_lin = 0.01*((20.0/0.01).^(x_lin/(41 - 1.0)));

x_lin = unique(x_lin);
y_lin = unique(y_lin);

% bad_id = (25 > avg_pwd);
% avg_pwd(bad_id) = NaN;
% avg_freq(bad_id) = NaN;
% avg_MFR_I(bad_id) = NaN;

% h = imagesc(avg_pwd);set(gca,'YDir','normal');set(h,'alphadata',~isnan(avg_pwd)) 

% x_lin_correct_unit = x_lin./wb_A_i.*(1e-9).*(1e+3);
% y_lin_correct_unit = y_lin./wb_A_i.*(1e-9).*(1e+3);
x_lin_correct_unit = x_lin;
y_lin_correct_unit = y_lin;

% my_imagesc_2D(x_lin, y_lin, avg_freq, ...
%     'is_XLog', 1, 'is_YLog', 0, ...
%     'CLimBegin', 20, 'CLimEnd', 100, ...
%     'is_ShowEdge', 0);

my_imagesc_2D(x_lin_correct_unit, y_lin_correct_unit, avg_pwd, ...
    'is_XLog', 0, 'is_YLog', 0, ...
    'CLimBegin', NaN, 'CLimEnd', NaN, ...
    'is_ShowEdge', 0);
colorbar

% set(gca, 'XTick', [1e-4 1e-3 1e-2 1e-1],  'XTickLabel',{'';'';'';''});
% set(gca, 'YTick', [0 0.005 0.01 0.015 0.02], 'YTickLabel',{'';'';'';'';''});
hTitle = title('Freq. [Hz] (\sigma_I = 0.5 mV & I_{0,I} = 1.1 \muA/cm^2)');
hXLabel = xlabel('g_{I->I} [mS/cm^2]');
hYLabel = ylabel('g_{GJ} [mS/cm^2]');
% axis tight
axis square
box on
% view(0,90);
% xlim([0.01 0.4]);
% ylim([0.001 0.01]);

make_me_pretty(gcf, ...
    gca, 15*1.5, ...
    hTitle, 30, ...
    hXLabel, 30, ...
    hYLabel, 30, ...
    [], 30, ...
    [], 12, ...
    [], 12)

maximize_a_fig(gcf);
% savefig('gII_vs_gGJII_sigmaWNI0_5_for_paper3', 'eps');
% 
% get_colorbar(0, 20, 100, 5)


end


