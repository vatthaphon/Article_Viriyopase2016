function my_plot
clc;clear all;close all;
format long

fig_id = 1;

Perc = 'X';

dir = 'E:\paper2_Raoul\Sim_network_of_other_people\data\nw_hh\EEEIIEII\CA1_IappEx_sigmaWNEx_pII0.3_pIIGJ0.004\v2\';
load(strcat(dir,strcat(Perc, '.mat')))

x_lin = unique(x_lin);m_xlabel = 'I_{0,E} [\muA/cm^2]';
y_lin = unique(y_lin);m_ylabel = '\sigma_{E} [mV]';

b_PowerFreq_dt = (E_PowerFreq_dt < 0.05) | isnan(E_PowerFreq_dt);
b_E_MFR_dt = (MFR_dt(:, 1) < 0) | (MFR_dt(:, 1) > 100);
b_I_MFR_dt = (MFR_dt(:, 2) < 0) | (MFR_dt(:, 2) > 100);
b_Freq_dt = (Freq_dt(:, 1) < 0) | (Freq_dt(:, 1) > 100);

b_cond = b_Freq_dt | b_PowerFreq_dt | b_E_MFR_dt | b_I_MFR_dt;

PowerFreq_dt(b_cond, 1) = NaN;
MFR_dt(b_cond, 1) = NaN;
MFR_dt(b_cond, 2) = NaN;
Freq_dt(b_cond, 1) = NaN;
E_Freq_dt(b_cond, 1) = NaN;
kappa_dt(b_cond, 1) = NaN;
kappa_dt(b_cond, 2) = NaN;
t_lags_dt(b_cond, 1) = NaN;

figure(fig_id);fig_id=fig_id+1;
% maximize_a_fig(gcf);
tmp = reshape(PowerFreq_dt(:, 1), N_y, N_x);m_title = 'Power';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
colorbar
% view(0,90);

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

%savefig('Power','eps');

figure(fig_id);fig_id=fig_id+1;
% maximize_a_fig(gcf);
subplot(2, 3, 1)
tmp = reshape(MFR_dt(:, 1), N_y, N_x);m_title = 'Mean firing rate of E-cells';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
% caxis([0.0 1.0]);
colorbar
%view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

subplot(2, 3, 2)
tmp = reshape(MFR_dt(:, 2), N_y, N_x);m_title = 'Mean firing rate of I-cells';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
caxis([0.0 50.0]);
colorbar
%view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

subplot(2, 3, 3)
tmp = reshape(Freq_dt(:, 1), N_y, N_x);m_title = 'Frequency';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
caxis([0.0 50.0]);
colorbar
% view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

subplot(2, 3, 4)
tmp = reshape(kappa_dt(:, 1), N_y, N_x);m_title = '\kappa of E-cells';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
% caxis([0.0 1.0]);
colorbar
%view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

subplot(2, 3, 5)
tmp = reshape(kappa_dt(:, 2), N_y, N_x);m_title = '\kappa of I-cells';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
% caxis([0.0 1.0]);
colorbar
%view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

subplot(2, 3, 6)
t_lags_dt((t_lags_dt<=-5.0) | (t_lags_dt>=50.0)) = NaN;
tmp = reshape(t_lags_dt(:, 1), N_y, N_x);m_title = 'Lag';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);

axis square
axis tight
% caxis([-1.0 5.0]);
colorbar
view(0,90);
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

%savefig('Analysis','eps');

% figure(fig_id);fig_id=fig_id+1;
% tmp = reshape(MFR_dt(:, 1)./MFR_dt(:, 2), N_y, N_x);m_title = 'Mean firing rate of E-cells';
% h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
% xlabel(m_xlabel);ylabel(m_ylabel);title(m_title);
% 
% axis square
% % caxis([0.0 1.0]);
% colorbar
% %view(0,90);
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)


nw_A_e = (21590*1e-8); % [cm^2]
x_lin = x_lin./nw_A_e.*(1e-12).*(1e+6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fig_id);fig_id=fig_id+1;
tmp = reshape(E_Freq_dt(:, 1), N_y, N_x);m_title = 'Freq. [Hz]';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
caxis([0.0 70.0]);
colorbar
% % % colorbar('YTick',[0 10 20 30 40], 'YTickLabel', {'';'';'';'';''})    
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0 5 10 15 20 25],'YTickLabel',{'';'';'';'';'';''});

axis square
axis tight

mXlabel = xlabel(m_xlabel);
mYlabel = ylabel(m_ylabel);
mTitle = title(m_title);

make_me_pretty(gcf, ...
    gca, 30, ...
    mTitle, 30, ...
    mXlabel, 30, ...
    mYlabel, 30, ...
    [], 12, ...
    [], 12, ...
    [], 12)
maximize_a_fig(gcf);
savefig('IappE_vs_SigmaWNE_Freq_HigherRange_for_paper3', 'eps');

figure(fig_id);fig_id=fig_id+1;
tmp = reshape(MFR_dt(:, 1), N_y, N_x);m_title = 'MFR_E [spks/s]';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
caxis([0.0 70.0]);
colorbar
% % colorbar('YTick',[0 10 20 30 40], 'YTickLabel', {'';'';'';'';''})    
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0 5 10 15 20 25],'YTickLabel',{'';'';'';'';'';''});

axis square
axis tight

mXlabel = xlabel(m_xlabel);
mYlabel = ylabel(m_ylabel);
mTitle = title(m_title);

make_me_pretty(gcf, ...
    gca, 30, ...
    mTitle, 30, ...
    mXlabel, 30, ...
    mYlabel, 30, ...
    [], 12, ...
    [], 12, ...
    [], 12)
maximize_a_fig(gcf);
savefig('IappE_vs_SigmaWNE_MFRe_HigherRange_for_paper3', 'eps');

figure(fig_id);fig_id=fig_id+1;
tmp = reshape(MFR_dt(:, 2), N_y, N_x);m_title = 'MFR_I [spks/s]';
h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
caxis([0.0 70.0]);
colorbar
% % colorbar('YTick',[0 10 20 30 40], 'YTickLabel', {'';'';'';'';''})    
% set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0 5 10 15 20 25],'YTickLabel',{'';'';'';'';'';''});

axis square
axis tight

mXlabel = xlabel(m_xlabel);
mYlabel = ylabel(m_ylabel);
mTitle = title(m_title);

make_me_pretty(gcf, ...
    gca, 30, ...
    mTitle, 30, ...
    mXlabel, 30, ...
    mYlabel, 30, ...
    [], 12, ...
    [], 12, ...
    [], 12)
maximize_a_fig(gcf);
savefig('IappE_vs_SigmaWNE_MFRi_HigherRange_for_paper3', 'eps');

% figure(fig_id);fig_id=fig_id+1;
% tmp = reshape(kappa_dt(:, 1), N_y, N_x);
% h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
% caxis([0.0 1.0]);
% colorbar
% % % colorbar('YTick',[0 0.25 0.5 0.75 1], 'YTickLabel', {'';'';'';'';''})    
% % set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% % set(gca,'YTick',[0 5 10 15 20],'YTickLabel',{'';'';'';'';''});
% 
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% savefig('IappE_vs_SigmaWNE_Ke', 'eps');
% 
% figure(fig_id);fig_id=fig_id+1;
% tmp = reshape(kappa_dt(:, 2), N_y, N_x);
% h = imagesc(x_lin, y_lin, tmp);set(gca,'YDir','normal');set(h,'alphadata',~isnan(tmp)) 
% caxis([0.0 1.0]);
% colorbar
% % colorbar('YTick',[0 0.25 0.5 0.75 1], 'YTickLabel', {'';'';'';'';''})    
% % set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',{'';'';'';'';''});
% % set(gca,'YTick',[0 5 10 15 20],'YTickLabel',{'';'';'';'';''});
% 
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% savefig('IappE_vs_SigmaWNE_Ki', 'eps');

% get_colorbar(1, 5)
end


