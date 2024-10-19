function showResultFromC_SaveRaster
clc;clear all;
close all;
format long

figure_id=1;

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi3_gIIGJx4_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi3_gIIGJx4_rndV0', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('C:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi3_gIIGJx3.6_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi3_gIIGJx3_6_rndV0', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi37_gIIGJx1.8_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi37_gIIGJx1_8_rndV0_39_28', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi3_gIIGJx0_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi3_gIIGJx0_rndV0_0_62', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi15_gIIGJx3.6_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi15_gIIGJx3_6_rndV0_58_58', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
% FN='NWB_II_gIIxi22_gIIGJx3.6_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
% FN_raster = strcat('gIIxi22_gIIGJx3_6_rndV0_76_38', '_Rasterogram');
% tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

% sub=strcat('E:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1', '\');
% FN='NWB_II_gIIxi0_gIIGJx0_0_0_1_lIIx0.6_sigmaWNI0.5_pIIGJ0_IappI200_rndV0_dt0.01_tEnd500_Ni1000';
% FN_raster = strcat('gIIxi0_gIIGJx0_rndV0', '_Rasterogram');
% tRasBegin=400;tRasEnd=500;E_color = 'r';I_color = 'b';

sub=strcat('C:\paper2_Raoul\Sim_network_of_other_people\data\hh\II\NWB_gIIxi_gIIGJx_0.01_20_41_sigmaWNI0.5_IappI200\v1\data', '\');
FN='NWB_II_gIIxi28_gIIGJx0.9_0.01_20_41_lIIx0.6_sigmaWNI0.5_pIIGJ0.004_IappI200_rndV0_dt0.01_tEnd2000_Ni1000';
FN_raster = strcat('gIIxi28_gIIGJx0_9_rndV0', '_Rasterogram');
tRasBegin=1900;tRasEnd=2000;E_color = 'r';I_color = 'b';

traj_ext='_Traj';
LE_ext='_LEs';  
LE_End_ext='_LEs_End';
LE_End_Pert_ext = '_LEs_End_Pert';
LE_Pert_ext = '_LEs_Pert';

N_LE = 2;tShowLEPert = 100;

LSpec_max_connected_I   = 'b';
LSpec_min_connected_I   = 'b--';
LSpec_max_connected_E   = 'r';
LSpec_min_connected_E   = 'r--';

isShowRastergram    = 1;isShowLabels = 1;
isShowPopCohMeasure = 0;
isShowMFR           = 1;isShowMFRFig = 0;
isComputeTimeLag    = 0;isShowComputeTimeLag=1;max_power_tol_percent_ComputeTimeLag = 0;
isShowFreqEachPop   = 0;
isShowProbabilityDischarge = 0;
isShowFreq          = 0;isShowFreqFig = 1;max_power_tol_percent = 50;
isShowSynIndex      = 0;
isShowVtraj         = 0;isShowVtrajZoom = 1;
isShowVtrajFigure   = 0;
isShowLE            = 0;
isShowLEEnd         = 0;
isSaveVtrajLE       = 0;
isCalEndPertLE      = 0;
isShowLEPert        = 0;

% Rasterogram
LineW = 1.5;
LineS = '-';
isShowAllOsc = 0;Max_N_i = 100;Max_N_e = 4*Max_N_i;N_showOSC=100;
vert_bar_size = 0.5*3.0*0.5*1*2.5;
isSaveRastergram = 0;
isShowPattSeq = 1;

% Population Coherence Measure
tau_popcohmeasure   = 1; % [] = ms.

% isShowProbabilityDischarge
isSaveProbabilityDischarge = 0;

% Lyapunov exponents
isLELogScale = 1;
isSaveLE = 0;
isShowEvolution = 0;
isShowDiscrete = 1;
r_end_LEs_levels = 1;
g_end_LEs_levels = 0;
b_end_LEs_levels = 0;
pert_size = 1e-6;

% isCalEndPertLE
i_LEs_pert = 120;
N_state_for_each_osc = 12;
v_FaceColor = 'r';
h_FaceColor = 'g';
n_FaceColor = 'b';
else_FaceColor = 'w';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID  =fopen(strcat(sub,FN), 'r');
if (0 <= fileID)
    Data    =fread(fileID, inf, 'double');fclose(fileID);
    
    [N_i N_e dt Nt SpkProfile p_EE p_EI p_IE p_II tau_m_on_E tau_m_on_I RmIe_E RmIe_I] = get_necessary_Data(Data);
    
    tSimBegin           = 100;
    tSimEnd             = dt*(Nt-1);
    
    resolution_dt       = 1; %[ms]
    resolution_dt_time_lag = dt; %[ms]
    resolution_dt_pdischarge = 0.1; %[ms]
    
    tMFRBegin           = 0;
    tMFREnd             = 50;    
    
    tFreqBegin          = tSimBegin;
    tFreqEnd            = tSimEnd;
    
    tTrajBegin          = tSimBegin;
    tTrajEnd            = tSimEnd;
    
    % tBegin_LEs          = tSimEnd*20/100;
    tBegin_LEs          = tSimEnd*0/100;
    tEnd_LEs            = tSimEnd;    
    
    sim_time=dt*(Nt-1);
    N_oscs=N_i+N_e;
else
    display(strcat('Cannot open:', strcat(sub,FN,traj_ext)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowPopCohMeasure == 1)
    [rows,cols] = find(SpkProfile == -1);
    
    eff_N_oscs = 100;
    kappa_ij = NaN(eff_N_oscs, eff_N_oscs);
    
    for i = 1:eff_N_oscs
        first_spk_times = [];
        
        if (i == 1)
            if (0 < (rows(i,1) - 1))
                first_spk_times = SpkProfile(1:(rows(i,1)-1),1);
                index=(tSimBegin <= first_spk_times) & (first_spk_times <= tSimEnd);
                first_spk_times = first_spk_times(index, 1);
            end
        else
            if ((rows(i-1,1) + 1) <= (rows(i,1) - 1))
                first_spk_times = SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1);
                index=(tSimBegin <= first_spk_times) & (first_spk_times <= tSimEnd);
                first_spk_times = first_spk_times(index, 1);
            end
        end
        
        for j = i:eff_N_oscs
            second_spk_times = [];
            
            if (j == 1)
                if (0 < (rows(j,1) - 1))
                    second_spk_times = SpkProfile(1:(rows(j,1)-1),1);
                    index=(tSimBegin <= second_spk_times) & (second_spk_times <= tSimEnd);
                    second_spk_times = second_spk_times(index, 1);
                end
            else
                if ((rows(j-1,1) + 1) <= (rows(j,1) - 1))
                    second_spk_times = SpkProfile((rows(j-1,1)+1):(rows(j,1)-1),1);
                    index=(tSimBegin <= second_spk_times) & (second_spk_times <= tSimEnd);
                    second_spk_times = second_spk_times(index, 1);
                end
            end
            
            [n_first_spk_times ,unused] = hist(first_spk_times,tSimBegin:tau_popcohmeasure:tSimEnd);
            [n_second_spk_times,unused] = hist(second_spk_times,tSimBegin:tau_popcohmeasure:tSimEnd);
            
            n_first_spk_times = (0 < n_first_spk_times).*1;
            n_second_spk_times = (0 < n_second_spk_times).*1;
            
            if ((max(n_first_spk_times) == 0) && (max(n_second_spk_times) == 0))
                kappa_ij(i,j) = 0.0;
                kappa_ij(j,i) = kappa_ij(i,j);
            else
                if (max(n_first_spk_times) == 0)
                    n_first_spk_times(1,1) = 1;
                end
                
                if (max(n_second_spk_times) == 0)
                    n_second_spk_times(1,1) = 1;
                end
                
                kappa_ij(i,j) = dot(n_first_spk_times, n_second_spk_times)/sqrt(sum(n_first_spk_times)*sum(n_second_spk_times));
                kappa_ij(j,i) = kappa_ij(i,j);
            end
        end
    end
    
    kappa = mean(mean(kappa_ij));
    
    display(strcat('Population Coherence measure:',num2str(kappa)));
end

if (isShowRastergram == 1)
    %% Plot rasterogram
    [rows,cols]=find(SpkProfile==-1);
    
    if (isShowAllOsc == 1)
        Show_ID=1:N_oscs;
    else
        if (N_oscs < 100)
            Show_ID=1:N_oscs;
        else
            if ((N_i == 0) || (N_e == 0))
                Show_ID=1:N_showOSC;
            else
                Show_ID=[1:1:Max_N_i (N_i + 1):1:(N_i + Max_N_e)];
            end
        end
    end
    
    running_id=1;
    figure(figure_id);hold on
    for i = Show_ID
        if (i<=N_i)
            if (i==1)
                if (0<(rows(i,1)-1))
                    tmp=SpkProfile(1:(rows(i,1)-1),1);
                    index=(tRasBegin<=tmp) & (tmp<=tRasEnd);
                    if (size(tmp(index,1),1)~=0)
                        n_x = size(tmp(index,1),1);
                        X = [tmp(index,1)'          ;tmp(index,1)'];
                        Y = [ones(1, n_x).*running_id     ;ones(1, n_x).*(running_id + vert_bar_size)];
                        line(X,Y,'Color',I_color,'LineWidth',LineW, 'LineStyle', LineS);
                    end
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    tmp=SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1);
                    index=(tRasBegin<=tmp) & (tmp<=tRasEnd);
                    if (size(tmp(index,1),1)~=0)
                        n_x = size(tmp(index,1),1);
                        X = [tmp(index,1)'          ;tmp(index,1)'];
                        Y = [ones(1, n_x).*running_id     ;ones(1, n_x).*(running_id + vert_bar_size)];
                        line(X,Y,'Color',I_color,'LineWidth',LineW, 'LineStyle', LineS);
                    end
                end
            end
        else
            if (i==1)
                if (0<(rows(i,1)-1))
                    tmp=SpkProfile(1:(rows(i,1)-1),1);
                    index=(tRasBegin<=tmp) & (tmp<=tRasEnd);
                    if (size(tmp(index,1),1)~=0)
                        n_x = size(tmp(index,1),1);
                        X = [tmp(index,1)'          ;tmp(index,1)'];
                        Y = [ones(1, n_x).*running_id     ;ones(1, n_x).*(running_id + vert_bar_size)];
                        line(X,Y,'Color',E_color,'LineWidth',LineW, 'LineStyle', LineS);
                    end
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    tmp=SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1);
                    index=(tRasBegin<=tmp) & (tmp<=tRasEnd);
                    if (size(tmp(index,1),1)~=0)
                        n_x = size(tmp(index,1),1);
                        X = [tmp(index,1)'          ;tmp(index,1)'];
                        Y = [ones(1, n_x).*running_id     ;ones(1, n_x).*(running_id + vert_bar_size)];
                        line(X,Y,'Color',E_color,'LineWidth',LineW, 'LineStyle', LineS);
                    end
                end
            end
        end
        
        running_id=running_id+1;
    end
    h=figure(figure_id);hold on;

%     hXLabel = xlabel('time [ms]');
%     hYLabel = ylabel('ID');
%     axis tight
%     axis square
    box on
        set(gca, 'XTick', [1900 1950 2000],  'XTickLabel',{'';'';''});
        set(gca, 'YTick', [1 50 100], 'YTickLabel',{'';'';''});
%         set(gca, 'YTick', [1 5 10], 'YTickLabel',{'';'';''});
    

    if (isShowLabels == 0)
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);        
    end
        
    if (min(Show_ID) < size(Show_ID, 2))
        ylim([min(Show_ID)-1 size(Show_ID, 2)+1]);
    end
    xlim([tRasBegin tRasEnd]);
    
    grid off;
%     ylim([0.5 11]);        
%     grey = [0.8,0.8,0.8];
%     
%     osc_id_xx = 1;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 2;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 3;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 4;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 5;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 6;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 7;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 8;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 9;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
%     osc_id_xx = 10;
%     plot([1900, 2000],[osc_id_xx, osc_id_xx], 'Color', grey);
    
    
    %     if (isShowPattSeq == 1)
    %         fileID  =fopen(strcat(sub,FN,'_PattSeq'), 'r');
    %         if (0 <= fileID)
    %             Data    =fread(fileID, inf, 'double');fclose(fileID);
    %             display(strcat('Patt seed:', num2str(Data(1,1))));
    %             display(strcat('N of Patts:', num2str(Data(2,1))));
    %
    %             pattseqs = Data(4:end, 1);
    %             t_pattseqs = pattseqs(1:2:(end - 1), 1);
    %             which_pattseqs = pattseqs(2:2:end, 1);
    %
    %             max_Show_ID = max(Show_ID);
    %             for i = 1:size(t_pattseqs, 1);
    %                 X = [t_pattseqs(i, 1)   ; t_pattseqs(i, 1)];
    %                 Y = [0                  ; (max_Show_ID + 1)];
    %
    %                 if (which_pattseqs(i, 1) == 1)
    %                     line(X,Y,'Color', 'b','LineWidth',LineW, 'LineStyle', '--');
    %                 elseif (which_pattseqs(i, 1) == 2)
    %                     line(X,Y,'Color', 'r','LineWidth',LineW, 'LineStyle', '--');
    %                 end
    %
    %
    %             end
    %         end
    %     end
    
%     if (isSaveRastergram == 1)
%         h=gcf;
%         set(h,'PaperOrientation','landscape');
%         set(h,'PaperUnits','normalized');
%         set(h,'PaperPosition', [0 0 1 1]);
%         print(gcf, '-dpdf', strcat(FN,'_raster.pdf'));
%         %         saveas(h, strcat(FN,'_raster.fig'), 'fig');
%     end
%     
%     make_me_pretty(gcf, gca, [], hXLabel, hYLabel, [], [], [], strcat(FN, '_Rasterogram.eps'));
%        
    make_me_pretty(gcf, ...
        gca, 12, ...
        [], 12, ...
        [], 12, ...
        [], 12, ...
        [], 12, ...
        [], 12, ...
        [], 12);
        
    maximize_a_fig(gcf);
    m_savefig(FN_raster, 'eps');
    
    grid off
    figure_id=figure_id+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowMFR == 1)
    display('Cal. MFR:');
    %% Calculate the mean firing rate
    MFR=SpkProfile(((tSimBegin <= SpkProfile) & (SpkProfile <= tSimEnd)) | (SpkProfile == -1));
    
    [rows,cols]=find(SpkProfile==-1);
    
    MFR(1,1)=(rows(1)-1)/sim_time;
    MFR(2:N_oscs,1)=(rows(2:N_oscs)-rows(1:(N_oscs-1))-1)/sim_time;
    
    MFR=MFR.*1000;  % Convert the [ms] to [s] so that we have [Hz] unit.
    MFR_E=MFR(N_i+1:N_oscs,1);
    MFR_I=MFR(1:N_i,1);
    
    display(strcat('- Avg. of MFR_E:', num2str(mean(mean(MFR_E))), ' Hz', ', ISI:', num2str(1000/mean(mean(MFR_E))), ' [ms]'));
    display(strcat('- Avg. of MFR_I:', num2str(mean(mean(MFR_I))), ' Hz', ', ISI:', num2str(1000/mean(mean(MFR_I))), ' [ms]'));
    
    [n_E,xout_E] = hist(MFR_E,tMFREnd - tMFRBegin);
    [n_I,xout_I] = hist(MFR_I,tMFREnd - tMFRBegin);
    
%     non_zeros=(n_E~=0);
%     n_E=n_E(non_zeros);
%     xout_E=xout_E(non_zeros);
%     
%     non_zeros=(n_I~=0);
%     n_I=n_I(non_zeros);
%     xout_I=xout_I(non_zeros);
    
    n_E=n_E/N_oscs;
    n_I=n_I/N_oscs;
    
    if (isShowMFRFig == 1)
        figure(figure_id)
        bar(xout_E,n_E,'r');
        hXLabel = xlabel('Firing rate [Hz]');
        hYLabel = ylabel('Fraction of cells');
        hTitle = title('Firing rate [Hz] of E-cells');
        axis tight
        xlim([tMFRBegin tMFREnd]);
        
        make_me_pretty(gcf, ...
            gca, 12, ...
            hTitle, 12, ...
            hXLabel, 12, ...
            hYLabel, 12, ...
            [], 12, ...
            [], 12, ...
            [], 12);
        
        figure_id=figure_id+1;
        
        figure(figure_id)
        bar(xout_I,n_I,'b');
        hXLabel = xlabel('Firing rate [Hz]');
        hYLabel = ylabel('Fraction of cells');
        hTitle = title('Firing rate [Hz] of I-cells');
        axis tight
        xlim([tMFRBegin tMFREnd]);
        
        make_me_pretty(gcf, ...
            gca, 12, ...
            hTitle, 12, ...
            hXLabel, 12, ...
            hYLabel, 12, ...
            [], 12, ...
            [], 12, ...
            [], 12);
        
        figure_id=figure_id+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isComputeTimeLag == 1)
    display('Cal. Time lag:');
    
    %% Calculate the population rhythmic frequency
    all_spike_times=SpkProfile((SpkProfile~=-1)); % [s]
    all_spike_times=all_spike_times(((tSimBegin<=all_spike_times) & (all_spike_times<=tSimEnd))); % [s]
    R_t=hist(all_spike_times,[tSimBegin:resolution_dt_time_lag:tSimEnd]); % Initialize the instantaneous firing rate with 'resolution_dt_time_lag' small bins.
    R_t=R_t/(N_oscs*resolution_dt_time_lag*1e-3);   % [Hz]
    
    N_x_t=size([tSimBegin:resolution_dt_time_lag:tSimEnd],2);
    t = [tSimBegin:resolution_dt_time_lag:tSimEnd];
    
    Fs      = 1000/resolution_dt_time_lag;                      % Sampling frequency [Hz]
    L       = N_x_t;                         % Length of signal
    
    NFFT    = 2^nextpow2(L); % Next power of 2 from length of y
    Y       = fft(R_t-mean(R_t),NFFT)/L;

    f       = Fs/2*linspace(0,1,NFFT/2+1);    
    power_Y = 2*abs(Y(1:NFFT/2+1));
    
    [max_power, I] = max(power_Y);
    tmp = f(1, (max_power - max_power*max_power_tol_percent_ComputeTimeLag/100) <= power_Y);
    T = 1000/tmp(1,1);
%     display(strcat('- Period for cal. time lag:', num2str(T)));
        
    %%
    [rows,cols]=find(SpkProfile==-1);
    I_spike_times=[];E_spike_times=[];
    for i=1:N_oscs
        if (i<=N_i)
            if (i==1)
                if (0<(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        else
            if (i==1)
                if (0<(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        end
    end
    
    I_spike_times=I_spike_times(((tSimBegin<=I_spike_times) & (I_spike_times<=tSimEnd)));
    E_spike_times=E_spike_times(((tSimBegin<=E_spike_times) & (E_spike_times<=tSimEnd)));
    
    [n_I_spk,t_I_spk]=hist(I_spike_times,[tSimBegin:resolution_dt_time_lag:tSimEnd]);
    [n_E_spk,t_E_spk]=hist(E_spike_times,[tSimBegin:resolution_dt_time_lag:tSimEnd]);
    
    n_I_spk=n_I_spk/(N_i*resolution_dt_time_lag*1e-3);   % [Hz]
    n_E_spk=n_E_spk/(N_e*resolution_dt_time_lag*1e-3);   % [Hz]
    
    [c,lags] = xcorr(n_I_spk, n_E_spk, 'coeff');
    t_lags_lin = lags*resolution_dt_time_lag;
    
    % return maximum value and its index
    [d,f] = max(c);
    
    t_lags1 = t_lags_lin(f);
    t_lags2 = -1*t_lags1/abs(t_lags1)*(T - abs(t_lags_lin(f)));
    
    t_lags = t_lags1;
    if abs(t_lags2) <= abs(t_lags1)
        t_lags = t_lags2;
    end
    
    if (0 <= t_lags)
        display(strcat('- E spike before I:', num2str(t_lags), ' ms.'));
    else
        display(strcat('- I spike before E:', num2str(abs(t_lags)), ' ms.'));
    end
    
    if (isShowComputeTimeLag == 1)
        figure(figure_id)
        plot(t_lags_lin, c);
        
        figure_id = figure_id +1;
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowFreqEachPop == 1)
    [rows,cols]=find(SpkProfile==-1);
    I_spike_times=[];E_spike_times=[];
    for i=1:N_oscs
        if (i<=N_i)
            if (i==1)
                if (0<(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        else
            if (i==1)
                if (0<(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        end
    end
    
    I_spike_times=I_spike_times(((tSimBegin<=I_spike_times) & (I_spike_times<=tSimEnd)));
    E_spike_times=E_spike_times(((tSimBegin<=E_spike_times) & (E_spike_times<=tSimEnd)));
    
    [n_I_spk,t_I_spk]=hist(I_spike_times,[tSimBegin:resolution_dt:tSimEnd]);
    [n_E_spk,t_E_spk]=hist(E_spike_times,[tSimBegin:resolution_dt:tSimEnd]);
    ylabel_txt = '#Spikes within 1 ms';
    
%     n_I_spk=n_I_spk/(resolution_dt*1e-3);   % [Hz]
%     n_E_spk=n_E_spk/(resolution_dt*1e-3);   % [Hz]
%     ylabel_txt = 'Population rate [Hz]';    
    
%     n_I_spk=n_I_spk/(N_i*resolution_dt*1e-3);   % [Hz]
%     n_E_spk=n_E_spk/(N_e*resolution_dt*1e-3);   % [Hz]
%     ylabel_txt = 'Population rate [Hz]';
    
    figure(figure_id)
    ref_figure_id = figure_id;
    rows = 3;cols = 1;id = 1;
    subplot(rows,cols,id)
    plot(t_I_spk,n_I_spk,'-b*',t_E_spk,n_E_spk,'-r*');
    legend('Interneuron','Pyramidal');
    xlabel('[ms]');ylabel(ylabel_txt);
    xlim([tFreqBegin tFreqEnd]);
    id = id + 1;
    
    subplot(rows,cols,id)
    plot(t_I_spk,n_I_spk,'-b*');
    legend('Interneuron');
    xlabel('[ms]');ylabel(ylabel_txt);
    xlim([tFreqBegin tFreqEnd]);
    id = id + 1;
    
    subplot(rows,cols,id)
    plot(t_E_spk,n_E_spk,'-r*');
    legend('Pyramidal');
    xlabel('[ms]');ylabel(ylabel_txt);
    xlim([tFreqBegin tFreqEnd]);
    id = id + 1;
    
    figure_id = figure_id +1;
    
    figure(figure_id)
    plot(t_I_spk,n_I_spk,'-b*',t_E_spk,n_E_spk,'-r*');
    hLegend = legend('Interneuron','Pyramidal', 'Location', 'NorthWest');
    hXLabel = xlabel('[ms]');
    hYLabel = ylabel(ylabel_txt);
    xlim([900 1000]);
%     ylim([0 350]);
    
    make_me_pretty(gcf, ...
        gca, 12, ...
        hTitle, 12, ...
        hXLabel, 12, ...
        hYLabel, 12, ...
        [], 12, ...
        [], 12, ...
        [], 12);
    
    figure_id = figure_id +1;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowProbabilityDischarge == 1)
    [rows,cols]=find(SpkProfile==-1);
    I_spike_times=[];E_spike_times=[];
    for i=1:N_oscs
        if (i<=N_i)
            if (i==1)
                if (0<(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    I_spike_times=[I_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        else
            if (i==1)
                if (0<(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile(1:(rows(i,1)-1),1)'];
                end
            else
                if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                    E_spike_times=[E_spike_times SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1)'];
                end
            end
        end
    end
    
    I_spike_times=I_spike_times(((tSimBegin<=I_spike_times) & (I_spike_times<=tSimEnd)));
    E_spike_times=E_spike_times(((tSimBegin<=E_spike_times) & (E_spike_times<=tSimEnd)));
    
    [n_I_spk,t_I_spk]=hist(I_spike_times,[tSimBegin:resolution_dt_pdischarge:tSimEnd]);
    [n_E_spk,t_E_spk]=hist(E_spike_times,[tSimBegin:resolution_dt_pdischarge:tSimEnd]);
    
    %% Determine probability of discharge by using the peak spiking times of the pyramidal cells of each blob as a reference
    n_t_E_spk = size(t_E_spk,2);
    cnt_E_blobs = 0;
    t_E_blobs_durations = [];
    t_E_blobs_begins = [];
    i_t_E_blobs_begins = [];
    t_E_blobs_ends = [];
    i_t_E_blobs_ends = [];
    abs_t_E_peak_each_blobs = [];
    rel_t_E_spk_times = [];
    i = 1;
    while (i <= n_t_E_spk)
        if (0 < n_E_spk(1,i))
            cnt_E_blobs = cnt_E_blobs + 1;
            t_E_blob_begin = t_E_spk(1,i);
            t_E_blobs_begins = [t_E_blobs_begins t_E_blob_begin];
            i_t_E_blobs_begins = [i_t_E_blobs_begins i];
            n_in_E_blob = [];
            t_in_E_blob = [];
            %% Get the blob
            cnt_con_zeros = 0;
            while ((cnt_con_zeros <= 50) && (i <= n_t_E_spk))
                if (n_E_spk(1,i) == 0)
                    cnt_con_zeros = cnt_con_zeros + 1;
                else
                    n_in_E_blob = [n_in_E_blob n_E_spk(1,i)];
                    t_in_E_blob = [t_in_E_blob t_E_spk(1,i)];
                    
                    t_E_blob_end = t_E_spk(1,i);
                    i_t_E_blob_end = i;
                    cnt_con_zeros = 0;
                end
                i = i + 1;
            end
            [C,I] = max(n_in_E_blob);
            t_max = t_in_E_blob(1, I);
            
            exp_spk_times = [];
            nn = size(n_in_E_blob,2);
            for ii=1:nn
                nnn = n_in_E_blob(1,ii);
                for iii=1:nnn
                    exp_spk_times = [exp_spk_times t_in_E_blob(1,ii)];
                end
            end
            rel_spk_times = exp_spk_times - t_max;
            
            
            abs_t_E_peak_each_blobs = [abs_t_E_peak_each_blobs t_max];
            rel_t_E_spk_times = [rel_t_E_spk_times rel_spk_times];
            
            t_E_blobs_ends = [t_E_blobs_ends t_E_blob_end];
            i_t_E_blobs_ends = [i_t_E_blobs_ends i_t_E_blob_end];
            
            t_E_blobs_durations = [t_E_blobs_durations (t_E_blob_end - t_E_blob_begin)];
        else
            i = i + 1;
        end
    end
    
    if (isShowFreqEachPop == 1)
        % Draw blob lines
        figure(ref_figure_id)
        rows = 3;cols = 1;id = 1;
        subplot(rows,cols,id)
        
        n_t_E_blobs_begins = size(t_E_blobs_begins, 2);
        x = [t_E_blobs_begins          ;t_E_blobs_begins];
        y = [ones(1, n_t_E_blobs_begins).*0     ;ones(1, n_t_E_blobs_begins).*40];
        line(x,y,'color','r');
        n_t_E_blobs_ends = size(t_E_blobs_ends, 2);
        x = [t_E_blobs_ends          ;t_E_blobs_ends];
        y = [ones(1, n_t_E_blobs_ends).*0     ;ones(1, n_t_E_blobs_ends).*40];
        line(x,y,'color','r');
        id = id + 1;
        
        subplot(rows,cols,id)
        id = id + 1;
        
        subplot(rows,cols,id)
        n_t_E_blobs_begins = size(t_E_blobs_begins, 2);
        x = [t_E_blobs_begins          ;t_E_blobs_begins];
        y = [ones(1, n_t_E_blobs_begins).*0     ;ones(1, n_t_E_blobs_begins).*40];
        line(x,y,'color','b');
        n_t_E_blobs_ends = size(t_E_blobs_ends, 2);
        x = [t_E_blobs_ends          ;t_E_blobs_ends];
        y = [ones(1, n_t_E_blobs_ends).*0     ;ones(1, n_t_E_blobs_ends).*40];
        line(x,y,'color','r');
        id = id + 1;
    end
    
    n_t_I_spk = size(t_I_spk,2);
    cnt_I_blobs = 0;
    t_I_blobs_durations = [];
    t_I_blobs_begins = [];
    i_t_I_blobs_begins = [];
    t_I_blobs_ends = [];
    i_t_I_blobs_ends = [];
    rel_t_I_spk_times = [];
    i = 1;
    while ((i <= n_t_I_spk) && (cnt_I_blobs < cnt_E_blobs))
        if (0 < n_I_spk(1,i))
            cnt_I_blobs = cnt_I_blobs + 1;
            t_I_blob_begin = t_I_spk(1,i);
            t_I_blobs_begins = [t_I_blobs_begins t_I_blob_begin];
            i_t_I_blobs_begins = [i_t_I_blobs_begins i];
            n_in_I_blob = [];
            t_in_I_blob = [];
            %% Get the blob
            cnt_con_zeros = 0;
            while ((cnt_con_zeros <= 50) && (i <= n_t_I_spk))
                if (n_I_spk(1,i) == 0)
                    cnt_con_zeros = cnt_con_zeros + 1;
                else
                    n_in_I_blob = [n_in_I_blob n_I_spk(1,i)];
                    t_in_I_blob = [t_in_I_blob t_I_spk(1,i)];
                    
                    t_I_blob_end = t_I_spk(1,i);
                    i_t_I_blob_end = i;
                    cnt_con_zeros = 0;
                end
                i = i + 1;
            end
            t_max = abs_t_E_peak_each_blobs(1, cnt_I_blobs);
            exp_spk_times = [];
            nn = size(n_in_I_blob,2);
            for ii=1:nn
                nnn = n_in_I_blob(1,ii);
                for iii=1:nnn
                    exp_spk_times = [exp_spk_times t_in_I_blob(1,ii)];
                end
            end
            rel_spk_times = exp_spk_times - t_max;
            
            rel_t_I_spk_times = [rel_t_I_spk_times rel_spk_times];
            
            t_I_blobs_ends = [t_I_blobs_ends t_I_blob_end];
            i_t_I_blobs_ends = [i_t_I_blobs_ends i_t_I_blob_end];
            
            t_I_blobs_durations = [t_I_blobs_durations (t_I_blob_end - t_I_blob_begin)];
        else
            i = i + 1;
        end
    end
    
    if (isShowFreqEachPop == 1)
        % Draw blob lines
        figure(ref_figure_id)
        rows = 3;cols = 1;id = 1;
        subplot(rows,cols,id)
        
        n_t_I_blobs_begins = size(t_I_blobs_begins, 2);
        x = [t_I_blobs_begins          ;t_I_blobs_begins];
        y = [ones(1, n_t_I_blobs_begins).*0     ;ones(1, n_t_I_blobs_begins).*80];
        line(x,y,'color','b');
        n_t_I_blobs_ends = size(t_I_blobs_ends, 2);
        x = [t_I_blobs_ends          ;t_I_blobs_ends];
        y = [ones(1, n_t_I_blobs_ends).*0     ;ones(1, n_t_I_blobs_ends).*80];
        line(x,y,'color','b');
        id = id + 1;
        
        subplot(rows,cols,id)
        n_t_I_blobs_begins = size(t_I_blobs_begins, 2);
        x = [t_I_blobs_begins          ;t_I_blobs_begins];
        y = [ones(1, n_t_I_blobs_begins).*0     ;ones(1, n_t_I_blobs_begins).*80];
        line(x,y,'color','b');
        n_t_I_blobs_ends = size(t_I_blobs_ends, 2);
        x = [t_I_blobs_ends          ;t_I_blobs_ends];
        y = [ones(1, n_t_I_blobs_ends).*0     ;ones(1, n_t_I_blobs_ends).*80];
        line(x,y,'color','r');
        id = id + 1;
        
        subplot(rows,cols,id)
        id = id + 1;
    end
    
    % Check if the blob devisions are correct
    if (cnt_I_blobs ~= cnt_E_blobs)
        display('Warning !!!: cnt_I_blobs ~= cnt_E_blobs');
    else
        %         if (isCorrectCheckBlobs(n_E_spk, n_I_spk, i_t_E_blobs_begins, i_t_I_blobs_begins, i_t_E_blobs_ends, i_t_I_blobs_ends) == 1)
        [n_rel_t_E_spk_times,t_rel_t_E_spk_times]=hist(rel_t_E_spk_times,[-10:resolution_dt_pdischarge:10]);
        [n_rel_t_I_spk_times,t_rel_t_I_spk_times]=hist(rel_t_I_spk_times,[-10:resolution_dt_pdischarge:10]);
        n_rel_t_E_spk_times=n_rel_t_E_spk_times/cnt_E_blobs;
        n_rel_t_I_spk_times=n_rel_t_I_spk_times/cnt_I_blobs;
        
        %             max_n_rel_t_E_spk_times = max(max(n_rel_t_E_spk_times));
        max_n_rel_t_E_spk_times = 1.0;
        
        figure(figure_id);hold on
        plot(t_rel_t_E_spk_times,n_rel_t_E_spk_times./max_n_rel_t_E_spk_times,'r');
        plot(t_rel_t_I_spk_times,n_rel_t_I_spk_times./max_n_rel_t_E_spk_times,'b');
        xlabel('time [ms]');
        ylabel('The number of spikes');
        %             legend('Distribution of spike times of pyramidal cells','Distribution of spike times of interneurons');
        box on
        
        if (isSaveProbabilityDischarge == 1)
            h=gcf;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            print(gcf, '-dpdf', strcat(FN,'_pdischarge.pdf'));
        end
        
        figure_id = figure_id +1;
        %         end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowFreq == 1)
    display('Cal. Freq. of Rhythm:');
    %% Calculate the population rhythmic frequency
    all_spike_times=SpkProfile((SpkProfile~=-1)); % [s]
    all_spike_times=all_spike_times(((tSimBegin<=all_spike_times) & (all_spike_times<=tSimEnd))); % [s]
    R_t=hist(all_spike_times,[tSimBegin:resolution_dt:tSimEnd]); % Initialize the instantaneous firing rate with 'resolution_dt' small bins.
    R_t=R_t/(N_oscs*resolution_dt*1e-3);   % [Hz]
    
    N_x_t=size([tSimBegin:resolution_dt:tSimEnd],2);
    t = [tSimBegin:resolution_dt:tSimEnd];
    
    Fs      = 1000/resolution_dt;                      % Sampling frequency [Hz]
    L       = N_x_t;                         % Length of signal
    
    NFFT    = 2^nextpow2(L); % Next power of 2 from length of y
    Y       = fft(R_t-mean(R_t),NFFT)/L;
    f       = Fs/2*linspace(0,1,NFFT/2+1);
        
    power_Y = 2*abs(Y(1:NFFT/2+1));
    [max_power, I] = max(power_Y);    
    tmp = f(1, (max_power - max_power*max_power_tol_percent/100) <= power_Y);
    display(strcat('- Range:',num2str((max_power - max_power*max_power_tol_percent/100)), '<='));
    display(strcat('- Freq.:',num2str(tmp(1,1)),' Hz (T =', num2str(1000/tmp(1,1)), ' ms.)'));
    
    display(strcat('- Stationary firing rate:',num2str(mean(mean(R_t))),' Hz'));
    
    if (isShowFreqFig == 1)
        figure(figure_id)
        row=3;col=1;id=1;
        
        subplot(row,col,id)
        plot(t, R_t,'-b*');
        xlabel('[ms]');ylabel('Population rate [Hz]');
        xlim([tFreqBegin tFreqEnd]);
        id=id+1;
        
        subplot(row,col,id)
        plot(f,2*abs(Y(1:NFFT/2+1)))    % Plot single-sided amplitude spectrum.
        title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
        id=id+1;
        
        subplot(row,col,id)
        plot(f,2*abs(Y(1:NFFT/2+1)))    % Plot single-sided amplitude spectrum.
        title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
        xlim([0 1000]);
        id=id+1;
        
        figure_id=figure_id+1;
        
        figure(figure_id)
        plot(f,2*abs(Y(1:NFFT/2+1)))    % Plot single-sided amplitude spectrum.
        hXLabel = xlabel('Frequency (Hz)');
        hYLabel = ylabel('|Y(f)|');
        xlim([0 200]);
        box on

        make_me_pretty(gcf, gca, [], hXLabel, hYLabel, [], [], [], strcat(FN, '_FreqSpec.eps'));
               
        
        figure_id=figure_id+1;                               
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowSynIndex == 1)
    %% Calculate MFR of the network
    MFR=SpkProfile(((tSimBegin <= SpkProfile) & (SpkProfile <= tSimEnd)) | (SpkProfile == -1));
    
    [rows,cols]=find(SpkProfile==-1);
    
    MFR(1,1)=(rows(1)-1)/sim_time;
    MFR(2:N_oscs,1)=(rows(2:N_oscs)-rows(1:(N_oscs-1))-1)/sim_time;
    
    MFR=MFR.*1000;  % Convert the [ms] to [s] so that we have [Hz] unit.
    MFR_E=MFR(N_i+1:N_oscs,1);
    MFR_I=MFR(1:N_i,1);
    
    mean_MFR = mean(mean([MFR_E;MFR_I]));
    
    %% Calculate the population rhythmic frequency
    all_spike_times=SpkProfile((SpkProfile~=-1)); % [s]
    all_spike_times=all_spike_times(((tSimBegin<=all_spike_times) & (all_spike_times<=tSimEnd)))/1000; % [s]
    R_t=hist(all_spike_times,[tSimBegin/1000:resolution_dt/1000:tSimEnd/1000]); % Initialize the instantaneous firing rate with 'resolution_dt' small bins.
    
    R_t=R_t/(N_oscs*resolution_dt*1e-3);   % [Hz]
    t = [tSimBegin/1000:resolution_dt/1000:tSimEnd/1000]; % [ms]
    
    figure(figure_id)
    subplot(2,1,1)
    plot(t,R_t,'b-*')
    [c1,lags1] = xcorr(R_t,R_t,'biased');
    subplot(2,1,2)
    plot(lags1, c1);
    c1 = max(c1);
    c1 = c1./(mean_MFR^2);
    display(strcat('Synchronization index:', num2str(c1)));
    
    figure_id=figure_id+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowVtraj == 1)
    fileID=fopen(strcat(sub,FN,traj_ext), 'r');
    
    if (0 <= fileID)
        N_data = 6;
        A = exist('ttt','var');
        if (A ~= 0)
            % Portion of trajectory
            fseek(fileID, 8*24*1/dt*ttt, 0);
            data=fread(fileID, 24*1/dt*ttt_length, 'double');fclose(fileID);
            
            ii=1;
            mul=N_data*4;
            V_max_connected_I=data((0:(1/dt*ttt_length-1)).*mul+ii,1);ii=ii+1;
            V_min_connected_I=data((0:(1/dt*ttt_length-1)).*mul+ii,1);ii=ii+1;
            V_max_connected_E=data((0:(1/dt*ttt_length-1)).*mul+ii,1);ii=ii+1;
            V_min_connected_E=data((0:(1/dt*ttt_length-1)).*mul+ii,1);ii=ii+1;
            
            t_lin = 0:(1/dt*ttt_length-1);
            
            if (isShowVtrajFigure == 1)
                figure(figure_id);hold on
                
                row=1;col=1;
                id=1;
                
                subplot(row,col,id);hold on
                plot(t_lin,V_max_connected_I,LSpec_max_connected_I,t_lin,V_min_connected_I,LSpec_min_connected_I,t_lin,V_max_connected_E,LSpec_max_connected_E,t_lin,V_min_connected_E,LSpec_min_connected_E);
                title('V');
                grid on;
                id=id+1;
                
                figure_id=figure_id+1;
            end
        else
            % Whole trajectory
            data=fread(fileID, inf, 'double');fclose(fileID);
            
             
            t_lin=linspace(0,dt*(Nt-1),Nt);
            t_lin_filtered_id = (tRasBegin <= t_lin) & (t_lin <= tRasEnd);
            t_lin = t_lin(t_lin_filtered_id);
            
            ii=1;
            mul=N_data*4;
            V_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;V_max_connected_I=V_max_connected_I(t_lin_filtered_id,1);
            V_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;V_min_connected_I=V_min_connected_I(t_lin_filtered_id,1);
            V_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;V_max_connected_E=V_max_connected_E(t_lin_filtered_id,1);
            V_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;V_min_connected_E=V_min_connected_E(t_lin_filtered_id,1);
            
            I_syn_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_max_connected_I=I_syn_max_connected_I(t_lin_filtered_id,1);
            I_syn_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_min_connected_I=I_syn_min_connected_I(t_lin_filtered_id,1);
            I_syn_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_max_connected_E=I_syn_max_connected_E(t_lin_filtered_id,1);
            I_syn_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_min_connected_E=I_syn_min_connected_E(t_lin_filtered_id,1);
            
            I_syn_from_I_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_max_connected_I=I_syn_from_I_max_connected_I(t_lin_filtered_id,1);
            I_syn_from_I_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_min_connected_I=I_syn_from_I_min_connected_I(t_lin_filtered_id,1);
            I_syn_from_I_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_max_connected_E=I_syn_from_I_max_connected_E(t_lin_filtered_id,1);
            I_syn_from_I_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_min_connected_E=I_syn_from_I_min_connected_E(t_lin_filtered_id,1);
            
            I_syn_from_E_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_max_connected_I=I_syn_from_E_max_connected_I(t_lin_filtered_id,1);
            I_syn_from_E_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_min_connected_I=I_syn_from_E_min_connected_I(t_lin_filtered_id,1);
            I_syn_from_E_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_max_connected_E=I_syn_from_E_max_connected_E(t_lin_filtered_id,1);
            I_syn_from_E_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_min_connected_E=I_syn_from_E_min_connected_E(t_lin_filtered_id,1);
            
            I_syn_from_E_ext_input_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_ext_input_max_connected_I=I_syn_from_E_ext_input_max_connected_I(t_lin_filtered_id,1);
            I_syn_from_E_ext_input_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_ext_input_min_connected_I=I_syn_from_E_ext_input_min_connected_I(t_lin_filtered_id,1);
            I_syn_from_E_ext_input_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_ext_input_max_connected_E=I_syn_from_E_ext_input_max_connected_E(t_lin_filtered_id,1);
            I_syn_from_E_ext_input_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_E_ext_input_min_connected_E=I_syn_from_E_ext_input_min_connected_E(t_lin_filtered_id,1);
            
            I_syn_from_I_ext_input_max_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_ext_input_max_connected_I=I_syn_from_I_ext_input_max_connected_I(t_lin_filtered_id,1);
            I_syn_from_I_ext_input_min_connected_I=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_ext_input_min_connected_I=I_syn_from_I_ext_input_min_connected_I(t_lin_filtered_id,1);
            I_syn_from_I_ext_input_max_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_ext_input_max_connected_E=I_syn_from_I_ext_input_max_connected_E(t_lin_filtered_id,1);
            I_syn_from_I_ext_input_min_connected_E=data((0:(Nt-1)).*mul+ii,1);ii=ii+1;I_syn_from_I_ext_input_min_connected_E=I_syn_from_I_ext_input_min_connected_E(t_lin_filtered_id,1);
            
            % Calculate statistic of current %
            if (abs(p_EE) < 1e-13)
                I_AMPA_E = 0.0;
            else
                I_AMPA_E = mean([mean(I_syn_from_E_max_connected_E), mean(I_syn_from_E_min_connected_E)]);
            end
            
            if (abs(p_EI) < 1e-13)
                I_AMPA_I = 0.0;
            else
                I_AMPA_I = mean([mean(I_syn_from_E_max_connected_I), mean(I_syn_from_E_min_connected_I)]);
            end
            
            if (abs(p_IE) < 1e-13)
                I_GABA_E = 0.0;
            else
                I_GABA_E = mean([mean(I_syn_from_I_max_connected_E), mean(I_syn_from_I_min_connected_E)]);
            end
            
            if (abs(p_II) < 1e-13)
                I_GABA_I = 0.0;
            else
                I_GABA_I = mean([mean(I_syn_from_I_max_connected_I), mean(I_syn_from_I_min_connected_I)]);
            end
            
            I_from_E_ext_syn_on_E = mean([mean(I_syn_from_E_ext_input_max_connected_E), mean(I_syn_from_E_ext_input_min_connected_E)]);
            I_from_E_ext_syn_on_I = mean([mean(I_syn_from_E_ext_input_max_connected_I), mean(I_syn_from_E_ext_input_min_connected_I)]);
            I_from_I_ext_syn_on_E = mean([mean(I_syn_from_I_ext_input_max_connected_E), mean(I_syn_from_I_ext_input_min_connected_E)]);
            I_from_I_ext_syn_on_I = mean([mean(I_syn_from_I_ext_input_max_connected_I), mean(I_syn_from_I_ext_input_min_connected_I)]);
            
            I_from_E_ext_on_E = mean([mean(I_syn_from_E_ext_input_max_connected_E), mean(I_syn_from_E_ext_input_min_connected_E)]) + RmIe_E;
            I_from_E_ext_on_I = mean([mean(I_syn_from_E_ext_input_max_connected_I), mean(I_syn_from_E_ext_input_min_connected_I)]) + RmIe_I;
            
            I_AMPA_tol_E = I_from_E_ext_on_E + I_AMPA_E;
            I_AMPA_tol_I = I_from_E_ext_on_I + I_AMPA_I;
            
            I_tol_E = RmIe_E + I_AMPA_E + I_GABA_E + I_from_E_ext_syn_on_E + I_from_I_ext_syn_on_E;
            I_tol_I = RmIe_I + I_AMPA_I + I_GABA_I + I_from_E_ext_syn_on_I + I_from_I_ext_syn_on_I;
            
            % Statistic information %
            %             display(strcat('I_AMPA_E/I_GABA_E:', num2str(I_AMPA_E/I_GABA_E)));
            %             display(strcat('I_AMPA_I/I_GABA_I:', num2str(I_AMPA_I/I_GABA_I)));
            %             display(strcat('I_AMPA_tol_E/I_GABA_E:', num2str(I_AMPA_tol_E/I_GABA_E)));
            %             display(strcat('I_AMPA_tol_I/I_GABA_I:', num2str(I_AMPA_tol_I/I_GABA_I)));
            display(strcat('I_tol_on_E:', num2str(I_tol_E)));
            display(strcat('- RmIe_E:', num2str(RmIe_E)));
            display(strcat('- I_AMPA_E:', num2str(I_AMPA_E)));
            display(strcat('- I_GABA_E:', num2str(I_GABA_E)));
            display(strcat('- I_from_E_ext_on_E:', num2str(I_from_E_ext_syn_on_E)));
            display(strcat('- I_from_I_ext_on_E:', num2str(I_from_I_ext_syn_on_E)));
            display(strcat('I_tol_on_I:', num2str(I_tol_I)));
            display(strcat('- RmIe_I:', num2str(RmIe_I)));
            display(strcat('- I_AMPA_I:', num2str(I_AMPA_I)));
            display(strcat('- I_GABA_I:', num2str(I_GABA_I)));
            display(strcat('- I_from_E_ext_on_I:', num2str(I_from_E_ext_syn_on_I)));
            display(strcat('- I_from_I_ext_on_I:', num2str(I_from_I_ext_syn_on_I)));
            
            if (isShowVtrajFigure == 1)
                figure(figure_id);hold on
                
                row=5;col=1;
                id=1;
                
                subplot(row,col,id);hold on
                plot(t_lin,V_max_connected_I,LSpec_max_connected_I,t_lin,V_min_connected_I,LSpec_min_connected_I,t_lin,V_max_connected_E,LSpec_max_connected_E,t_lin,V_min_connected_E,LSpec_min_connected_E);
                title('V');
                xlim([tRasBegin tRasEnd]);
                grid on;
                id=id+1;
                
                subplot(row,col,id);hold on
                plot(t_lin,I_syn_from_I_max_connected_I,LSpec_max_connected_I,t_lin,I_syn_from_I_min_connected_I,LSpec_min_connected_I,t_lin,I_syn_from_I_max_connected_E,LSpec_max_connected_E,t_lin,I_syn_from_I_min_connected_E,LSpec_min_connected_E);
                title('g_{syn\_from\_I}');
                xlim([tRasBegin tRasEnd]);
                grid on;
                id=id+1;
                
                subplot(row,col,id);hold on
                plot(t_lin,I_syn_from_E_max_connected_I,LSpec_max_connected_I,t_lin,I_syn_from_E_min_connected_I,LSpec_min_connected_I,t_lin,I_syn_from_E_max_connected_E,LSpec_max_connected_E,t_lin,I_syn_from_E_min_connected_E,LSpec_min_connected_E);
                title('g_{syn\_from\_E}');
                xlim([tRasBegin tRasEnd]);
                grid on;
                id=id+1;
                
                subplot(row,col,id);hold on
                plot(t_lin,I_syn_from_E_ext_input_max_connected_I,LSpec_max_connected_I, ...
                    t_lin,I_syn_from_E_ext_input_min_connected_I,LSpec_min_connected_I, ...
                    t_lin,I_syn_from_E_ext_input_max_connected_E,LSpec_max_connected_E, ...
                    t_lin,I_syn_from_E_ext_input_min_connected_E,LSpec_min_connected_E);
                
                title('g_{syn\_E\_ext\_input}');
                xlim([tRasBegin tRasEnd]);
                grid on;
                id=id+1;
                
                subplot(row,col,id);hold on
                plot(t_lin,I_syn_from_I_ext_input_max_connected_I,LSpec_max_connected_I, ...
                    t_lin,I_syn_from_I_ext_input_min_connected_I,LSpec_min_connected_I, ...
                    t_lin,I_syn_from_I_ext_input_max_connected_E,LSpec_max_connected_E, ...
                    t_lin,I_syn_from_I_ext_input_min_connected_E,LSpec_min_connected_E);
                
                title('g_{syn\_I\_ext\_input}');
                xlim([tRasBegin tRasEnd]);
                grid on;
                id=id+1;
                
                figure_id=figure_id+1;                
            end
            
            if (isShowVtrajZoom == 1)
                if (0 < N_i)
                    figure(figure_id);hold on
                    
                    plot(t_lin,V_max_connected_I,'-b');                    
                    grid on 
                    
                    ylim([-80 60]);
                    
                    hXLabel = xlabel('time [ms]');
                    hYLabel = ylabel('Voltage [mV]');
                    box on
                    
                    make_me_pretty(gcf, gca, [], hXLabel, hYLabel, [], [], [], strcat(FN, '_Vi.eps'));
                    
                    figure_id=figure_id+1;                    
                end
                                
                if (0 < N_e)
                    figure(figure_id);hold on
                    
                    plot(t_lin,V_max_connected_E,'-r');                    
                    grid on                                        
                    
                    ylim([-80 60]);
                    
                    hXLabel = xlabel('time [ms]');
                    hYLabel = ylabel('Voltage [mV]');
                    box on
                    
                    make_me_pretty(gcf, gca, [], hXLabel, hYLabel, [], [], [], strcat(FN, '_Ve.eps'));                    
                    
                    figure_id=figure_id+1;                                    
                end
                
                
            end
        end
    else
        display(strcat('Cannot open:', strcat(sub,FN,traj_ext)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isShowLEEnd == 1)
    fileID = fopen(strcat(sub, FN, LE_End_ext), 'r');
    if (0 <= fileID)
        %% Open the LE_End file
        data = fread(fileID, inf, 'double');fclose(fileID);
        N_LE = data(1, 1);
        data = data(2:end, 1);
        LEs_End = reshape(data, N_LE, []);
        
        for i = 1:N_LE
            display(strcat('LE=', num2str(LEs_End(1,i)), ' Avg LE=', num2str(LEs_End(2,i)), ' Avg Avg LE=', num2str(LEs_End(3,i))));
        end
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_End_ext)));
    end
end

if (isShowLE == 1)
    fileID = fopen(strcat(sub, FN, LE_ext), 'r');
    if (0 <= fileID)
        %% Open the LE file
        data = fread(fileID, inf, 'double');fclose(fileID);
        
        %% Get the number of LE and all LEs
        N_LE = data(1, 1);
        data = data(2:end, 1);
        LEs = reshape(data, N_LE, []);
        
        %% Record the Lyapunov exponents
        end_LEs = NaN(N_LE, 1);
        
        %% Get mask
        max_LEs = LEs(1,:);
        n_max_LEs = size(max_LEs, 2);
        
        filter_t = linspace(0, tEnd_LEs, n_max_LEs)';
        t_mask = (filter_t >= tBegin_LEs);
        
        %         %% Get the maxmimum LE
        %         max_LEs = LEs(1,:);
        %         n_max_LEs = size(max_LEs, 2);
        %
        %         filter_t = linspace(0, tEnd_LEs, n_max_LEs)';
        %         mask = (filter_t >= tBegin_LEs);
        %
        %         max_LEs = max_LEs(1, mask);
        %
        %         figure(200);hold on;
        %         plot(max_LEs)
        %
        %         % Determine evolution of avg. value of LE
        %         evo_max_LEs = NaN(1, size(max_LEs, 2));
        %         evo_max_LEs(1, 1) = max_LEs(1, 1);
        %         for i = 2:size(max_LEs, 2)
        %             evo_max_LEs(1, i) = (evo_max_LEs(1, i - 1)*(i - 1) + max_LEs(1, i))/i;
        %         end
        %
        %         % Decompose into positive and negative avg LEs
        %         posi_max_LEs = max_LEs(0 < max_LEs);
        %         nega_max_LEs = max_LEs(max_LEs < 0);
        %
        %         display(strcat('Avg. positive contribution to LEs:', num2str(mean(posi_max_LEs))));
        %         display(strcat('Avg. negative contribution to LEs:', num2str(mean(nega_max_LEs))));
        %         display(strcat('Mean of LEs:', num2str(mean(max_LEs))));
        %
        %         %% Display
        %         figure(figure_id);
        % %         row = 2; col = 1;
        %         row = 1; col = 1;
        %         id = 1;
        %
        %         n_max_LEs = size(max_LEs, 2);
        % %         subplot(row, col, id);hold on
        % %         plot(linspace(tBegin_LEs, tEnd_LEs, n_max_LEs), max_LEs, '-b');
        % % %         save('no_re11.mat','max_LEs');
        % %         if (isLELogScale == 1)
        % %             set(gca, 'xscale', 'log');
        % %         end
        % %         id = id + 1;
        %
        %         subplot(row, col, id);hold on
        %         plot(linspace(tBegin_LEs - tBegin_LEs, tEnd_LEs - tBegin_LEs, n_max_LEs), evo_max_LEs, '-b', 'LineWidth', 2);
        % %         save('re13.mat','evo_max_LEs');
        %         if (isLELogScale == 1)
        %             set(gca, 'xscale', 'log');
        %         end
        %         id = id + 1;
        %
        %         %% Determine the max LE by using average of the last portion
        %         end_max_LE = get_end_max_LE(evo_max_LEs);
        %         grid on
        %         display(strcat('Max. LE:', num2str(end_max_LE)));
        %         display(strcat('Last max. LE:', num2str(evo_max_LEs(1, end))));
        %         box on
        %         axis tight
        %         title(strcat('\lambda=',num2str(end_max_LE)));
        %         ylabel('The maximum Lyapunov exponent');
        %         xlabel('Relative time [ms]');
        %
        %         end_LEs(1, 1) = end_max_LE;
        
        %% Show the next LEs
        r_levels = linspace(1, 0, N_LE);
        b_levels = linspace(0, 1, N_LE);
        for i = 1:N_LE
            LE = LEs(i,:);
            LE = LE(1, t_mask);
            
            if (isShowDiscrete == 1)
                figure(figure_id);hold on
                plot(linspace(tBegin_LEs - tBegin_LEs, tEnd_LEs - tBegin_LEs, size(LE,2)), LE, 'b*');
                ylabel('The Local Lyapunov exponents');
                xlabel('Relative time [ms]');
                grid on;
                box on;
                axis tight;
            end
            
            if (isShowEvolution == 1)
                evo_LE = NaN(1, size(LE, 2));
                evo_LE(1, 1) = LE(1, 1);
                for j = 2:size(LE, 2)
                    evo_LE(1, j) = (evo_LE(1, j - 1)*(j - 1) + LE(1, j))/j;
                end
                
                figure(figure_id);hold on
                plot(linspace(tBegin_LEs - tBegin_LEs, tEnd_LEs - tBegin_LEs, size(evo_LE,2)), evo_LE, 'Color', [r_levels(1,i) 0 b_levels(1,i)], 'LineWidth', 2);
                if (isLELogScale == 1)
                    set(gca, 'xscale', 'log');
                end
                ylabel('The Lyapunov exponents');
                xlabel('Relative time [ms]');
                grid on;
                box on;
                axis tight;
            else
                evo_LE = LE;
            end
            
            end_LEs(i, 1) = get_end_max_LE(evo_LE)
        end
        
        if (isSaveLE == 1)
            h=gcf;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            print(gcf, '-dpdf', strcat(FN,'_LE.pdf'));
        end
        
        %         figure_id=figure_id+1;
        %
        %         figure(figure_id);hold on;
        % %         plot(sort(end_LEs, 'descend'), '*');
        %         plot(end_LEs, '*', 'Color', [r_end_LEs_levels g_end_LEs_levels b_end_LEs_levels]);
        %         ylabel('The Lyapunov exponents');
        %         xlabel('n/N');
        %         grid on;
        %         box on;
        figure_id=figure_id+1;
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_ext)));
    end
end

if (isSaveVtrajLE == 1)
    if (exist('V_max_connected_I','var') == 1)
        save(strcat(FN, '.mat'), 'dt', 'Nt', 'N_i', 'N_e', 'V_max_connected_I', 'V_min_connected_I', 'V_max_connected_E', 'V_min_connected_E', 'tau_m_on_E', 'tau_m_on_I');
    else
        display('V_max_connected_I doesnot exist');
    end
end

if (isCalEndPertLE == 1)
    fileID = fopen(strcat(sub, FN, LE_End_Pert_ext), 'r');
    if (0 <= fileID)
        %% Open the LE file
        data = fread(fileID, inf, 'double');fclose(fileID);
        
        %% Get the number of LE and all LEs
        N_LE = data(1, 1);
        N_state = data(2, 1);
        data = data(3:end, 1);
        LEs_pert = reshape(data, N_state, []);
        
        %         N_LE_selected = 1;
        %         r_levels = linspace(1, 0, N_LE_selected);
        %         b_levels = linspace(0, 1, N_LE_selected);
        %         for i = 1:1:N_LE_selected
        %             figure(figure_id);hold on;
        %             plot(abs(LEs_pert(:, i))./pert_size, '-*', 'Color', [r_levels(1,i) 0 b_levels(1,i)], 'LineWidth', 1);
        %         end
        
        x_pie = abs(LEs_pert(:, i_LEs_pert)')./pert_size;
        
        explode_pie = zeros(1, size(x_pie, 2));
        explode_pie(1, 1:N_state_for_each_osc:N_state) = 1;
        
        figure(figure_id);hold on;
        h = pie(x_pie,explode_pie);
        for i = 1:1:N_state
            if (rem(i, 12) == 1)
                set(h(2*(i - 1) + 1), 'FaceColor', v_FaceColor);
            elseif (rem(i, 12) == 2)
                set(h(2*(i - 1) + 1), 'FaceColor', h_FaceColor);
            elseif (rem(i, 12) == 3)
                set(h(2*(i - 1) + 1), 'FaceColor', n_FaceColor);
            else
                set(h(2*(i - 1) + 1), 'FaceColor', else_FaceColor);
            end
        end
        
        %% Axis for all-to-all connections (except to themself)
        %         x_osc_tick = [];
        %         for i = 10:10:100
        %             x_osc_tick = [x_osc_tick i*102 + 1];
        %         end
        
        %         axis tight
        %         box on
        %         grid on
        %         ylim([0 1])
        %         xlabel('State variables');
        %         ylabel('Normalized absolute value of perturbation');
        
        %         legend('Corresponding to the maximum Lyapunov exponent', 'Corresponding to the next Lyapunov exponent');
        
        figure_id=figure_id+1;
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_End_Pert_ext)));
    end
end

if (isShowLEPert == 1)
    fileID = fopen(strcat(sub, FN, LE_Pert_ext), 'r');
    
    %% Open the LE file.
    if (0 <= fileID)
        fseek(fileID, 2*8, 'bof');
        
        angle= NaN(Nt, 1);
        
        i = 1;
        while (feof(fileID) == 0)
            data = fread(fileID, N_LE*N_osc, 'double');
            
            if (size(data, 1) ~= 0)
                LE_perts = reshape(data, N_osc, []);
                
                LE_pert1 = LE_perts(:, 1);
                LE_pert2 = LE_perts(:, 2);
                
                angle(i, 1) = acos(dot(LE_pert1, LE_pert2)/(norm(LE_pert1)*norm(LE_pert2))).*180./pi;
                
                i = i + 1;
                
            end
        end
        
        fclose(fileID);
        
        angle = angle(isnan(angle) == 0);
        
        figure(figure_id);hold on;
        plot(linspace(0, Nt*dt, size(angle, 1)), angle);
        ylabel('Degree');
        xlabel('Time [ms]');
        
        figure_id=figure_id+1;
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_Pert_ext)));
    end
    
    fileID = fopen(strcat(sub, FN, LE_Pert_ext), 'r');
    %% Open the LE file and read the first tShowLEPert ms.
    if (0 <= fileID)
        data = fread(fileID, 2 + tShowLEPert*N_LE*(N_i + N_e), 'double');fclose(fileID);
        N_LE = data(1, 1);
        N_osc = data(2, 1);
        data = data(3:end, 1);
        
        LE_perts = reshape(data, N_osc*N_LE, []);
        
        no_of_LE_perts = size(LE_perts, 2);
        
        angle = NaN(no_of_LE_perts, 1);
        norm_diff = NaN(no_of_LE_perts, 1);
        
        for i = 1:no_of_LE_perts
            LE_pert = reshape(LE_perts(:,i), N_osc, []);
            
            angle(i, 1) = acos(dot(LE_pert(:, 1), LE_pert(:, 2))/(norm(LE_pert(:, 1))*norm(LE_pert(:, 2)))).*180./pi;
            
            %             norm_diff(i, 1) = norm(LE_pert(:, 1) - LE_pert(:, 2))/PERT_SIZE_LEs;
        end
        
        figure(figure_id)
        %         rows = 2; cols = 1; id = 1;
        
        %         subplot(rows, cols, id)
        plot(angle);id = id + 1;
        
        %         subplot(rows, cols, id)
        %         plot(norm_diff);id = id + 1;
        
        figure_id=figure_id+1;
        
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_Pert_ext)));
    end
    
    fileID = fopen(strcat(sub, FN, LE_Pert_ext), 'r');
    %% Open the LE file and read the last tShowLEPert ms.
    if (0 <= fileID)
        fseek(fileID, -tShowLEPert*N_LE*(N_i + N_e)*8, 'eof');
        data = fread(fileID, tShowLEPert*N_LE*(N_i + N_e), 'double');fclose(fileID);
        
        LE_perts = reshape(data, N_osc*N_LE, []);
        
        no_of_LE_perts = size(LE_perts, 2);
        
        angle = NaN(no_of_LE_perts, 1);
        norm_diff = NaN(no_of_LE_perts, 1);
        
        for i = 1:no_of_LE_perts
            LE_pert = reshape(LE_perts(:,i), N_osc, []);
            
            angle(i, 1) = acos(dot(LE_pert(:, 1), LE_pert(:, 2))/(norm(LE_pert(:, 1))*norm(LE_pert(:, 2)))).*180./pi;
            
            %             norm_diff(i, 1) = norm(LE_pert(:, 1) - LE_pert(:, 2))/PERT_SIZE_LEs;
        end
        
        figure(figure_id)
        %         rows = 2; cols = 1; id = 1;
        
        %         subplot(rows, cols, id)
        plot(angle);id = id + 1;
        
        %         subplot(rows, cols, id)
        %         plot(norm_diff);id = id + 1;
        
        figure_id=figure_id+1;
        
    else
        display(strcat('Cannot open:', strcat(sub, FN, LE_Pert_ext)));
    end
end

end

function val = get_end_max_LE(Evo_LE)
%     [x,resnorm] = lsqcurvefit(@myfun_exp_de, [1; 1], 1:size(Evo_LE, 2), Evo_LE)
%
%     val = myfun_exp_de(x, size(Evo_LE, 2));

begin = floor(size(Evo_LE, 2)*50/100);
val = mean(Evo_LE(1, (end - begin):end));
end

function F = myfun_exp_de(x, xdata)
F = x(2)./x(1).*exp(-xdata./x(1));
end

function [x_t,R_t]=cal_R_t(SpkProfile,N_oscs,tBegin_x_t,tEnd_x_t,N_x_t)
t=linspace(tBegin_x_t,tEnd_x_t,N_x_t)';
dt=(tEnd_x_t-tBegin_x_t)/(N_x_t-1);
R_t=zeros(N_x_t,1);
x_t=NaN(N_x_t,1);

[rows,cols]=find(SpkProfile==-1);
for j=1:N_x_t
    for i=1:N_oscs
        if (i==1)
            if (0<(rows(i,1)-1))
                tmp=SpkProfile(1:(rows(i,1)-1),1);
                index=(t(j,1)<=tmp) & (tmp<t(j,1)+dt);
                R_t(j,1)=R_t(j,1)+sum(index);
            end
        else
            if ((rows(i-1,1)+1)<=(rows(i,1)-1))
                tmp=SpkProfile((rows(i-1,1)+1):(rows(i,1)-1),1);
                index=(t(j,1)<=tmp) & (tmp<t(j,1)+dt);
                R_t(j,1)=R_t(j,1)+sum(index);
            end
        end
    end
    x_t(j,1)=t(j,1);
end

R_t=R_t/(dt*N_oscs)*1000;
end

function val=cmp(x, y, tol_eq)
% CMP Two-value comparison
%   val = cmp(x, y, tol_eq)
% Input
%   x           the first number.
%   y           the second number.
%   tol_eq      if the first and second numbers are different less than
%               tol_eq, we say that the two numbers are equal.
% Output
%   val         0   : two numbers are the same.
%               -1  : the first number is less than the second number.
%               1   : the first number is greater than the second number.

if  (abs(x-y)<tol_eq)
    val=0;
elseif (x<y)
    val=-1;
else
    val=1;
end

end


function [N_i N_e dt Nt SpkProfile p_EE p_EI p_IE p_II tau_m_on_E tau_m_on_I RmIe_E RmIe_I] = get_necessary_Data(Data)
% %%
% cnt=1;
% 
% g_syn_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% g_syn_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% g_syn_GABA_on_E=Data(cnt,1);cnt=cnt+2;
% g_syn_GABA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% v_syn_AMPA=Data(cnt,1);cnt=cnt+2;
% v_syn_GABA=Data(cnt,1);cnt=cnt+2;
% 
% tau_m_on_E=Data(cnt,1);cnt=cnt+2;
% tau_m_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_l_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_l_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% tau_l_GABA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_l_GABA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_r_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_r_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% tau_r_GABA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_r_GABA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_d_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_d_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% tau_d_GABA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_d_GABA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% V_rest=Data(cnt,1);cnt=cnt+2;
% V_th=Data(cnt,1);cnt=cnt+2;
% V_reset=Data(cnt,1);cnt=cnt+2;
% V_refrac=Data(cnt,1);cnt=cnt+2;
% tau_m=Data(cnt,1);cnt=cnt+2;
% C_m=Data(cnt,1);cnt=cnt+2;
% 
% V_rest=Data(cnt,1);cnt=cnt+2;
% V_th=Data(cnt,1);cnt=cnt+2;
% V_reset=Data(cnt,1);cnt=cnt+2;
% V_refrac=Data(cnt,1);cnt=cnt+2;
% tau_m=Data(cnt,1);cnt=cnt+2;
% C_m=Data(cnt,1);cnt=cnt+2;
% 
% N_i=Data(cnt,1);cnt=cnt+2;
% N_e=Data(cnt,1);cnt=cnt+2;
% N_ext=Data(cnt,1);cnt=cnt+2;
% 
% p_EE=Data(cnt,1);cnt=cnt+2;
% p_EI=Data(cnt,1);cnt=cnt+2;
% p_IE=Data(cnt,1);cnt=cnt+2;
% p_II=Data(cnt,1);cnt=cnt+2;
% 
% lambda_on_E=Data(cnt,1);
% cnt=cnt+2;
% lambda_on_I=Data(cnt,1);
% cnt=cnt+2;
% 
% v_syn_AMPA=Data(cnt,1);cnt=cnt+2;
% 
% g_ext_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% g_ext_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_m_on_I=Data(cnt,1);cnt=cnt+2;
% tau_m_on_E=Data(cnt,1);cnt=cnt+2;
% 
% tau_l_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_l_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_r_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_r_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% tau_d_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
% tau_d_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
% 
% RmIe_E=Data(cnt,1);cnt=cnt+2;
% RmIe_I=Data(cnt,1);cnt=cnt+2;
% 
% dt=Data(cnt,1);cnt=cnt+2;
% Nt=Data(cnt,1);cnt=cnt+2;
% seed=Data(cnt,1);cnt=cnt+2;
% isUsingSeed=Data(cnt,1);cnt=cnt+2;
% 
% SpkProfile=Data(cnt:end,1);
% 
% if (isUsingSeed ~= 1)
    cnt=1;
    
    g_syn_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    g_syn_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    g_syn_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    g_syn_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    v_syn_AMPA=Data(cnt,1);cnt=cnt+2;
    v_syn_GABA=Data(cnt,1);cnt=cnt+2;
    
    tau_m_on_E=Data(cnt,1);cnt=cnt+2;
    tau_m_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_l_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_l_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_l_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_l_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_r_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_r_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_r_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_r_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_d_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_d_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_d_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_d_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    V_rest=Data(cnt,1);cnt=cnt+2;
    V_th=Data(cnt,1);cnt=cnt+2;
    V_reset=Data(cnt,1);cnt=cnt+2;
    V_refrac=Data(cnt,1);cnt=cnt+2;
    tau_m=Data(cnt,1);cnt=cnt+2;
    C_m=Data(cnt,1);cnt=cnt+2;
    
    V_rest=Data(cnt,1);cnt=cnt+2;
    V_th=Data(cnt,1);cnt=cnt+2;
    V_reset=Data(cnt,1);cnt=cnt+2;
    V_refrac=Data(cnt,1);cnt=cnt+2;
    tau_m=Data(cnt,1);cnt=cnt+2;
    C_m=Data(cnt,1);cnt=cnt+2;
    
    N_i=Data(cnt,1);cnt=cnt+2;
    N_e=Data(cnt,1);cnt=cnt+2;
    N_ext=Data(cnt,1);cnt=cnt+2;
    
    p_EE=Data(cnt,1);cnt=cnt+2;
    p_EI=Data(cnt,1);cnt=cnt+2;
    p_IE=Data(cnt,1);cnt=cnt+2;
    p_II=Data(cnt,1);cnt=cnt+2;
    
    lambda_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    lambda_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    lambda_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    lambda_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    v_syn_ext_AMPA=Data(cnt,1);cnt=cnt+2;
    v_syn_ext_GABA=Data(cnt,1);cnt=cnt+2;
    
    g_ext_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    g_ext_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    g_ext_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    g_ext_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_m_ext_on_I=Data(cnt,1);cnt=cnt+2;
    tau_m_ext_on_E=Data(cnt,1);cnt=cnt+2;
    
    tau_l_ext_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_l_ext_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_l_ext_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_l_ext_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_r_ext_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_r_ext_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_r_ext_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_r_ext_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    tau_d_ext_AMPA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_d_ext_AMPA_on_I=Data(cnt,1);cnt=cnt+2;
    tau_d_ext_GABA_on_E=Data(cnt,1);cnt=cnt+2;
    tau_d_ext_GABA_on_I=Data(cnt,1);cnt=cnt+2;
    
    RmIe_E=Data(cnt,1);cnt=cnt+2;
    RmIe_I=Data(cnt,1);cnt=cnt+2;
    
    dt=Data(cnt,1);cnt=cnt+2;
    Nt=Data(cnt,1);cnt=cnt+2;
    seed=Data(cnt,1);cnt=cnt+2;
    isUsingSeed=Data(cnt,1);cnt=cnt+2;
    
    SpkProfile=Data(cnt:end,1);
% end

%% Display
display(strcat('N_e:', num2str(N_e)));
display(strcat('N_i:', num2str(N_i)));

display(strcat('RmIe_E:', num2str(RmIe_E)));
display(strcat('RmIe_I:', num2str(RmIe_I)));
display(strcat('Simulation time:', num2str(dt*Nt)));
display(strcat('dt:', num2str(dt)));
display(strcat('pEE:', num2str(p_EE)));
display(strcat('pEI:', num2str(p_EI)));
display(strcat('pIE:', num2str(p_IE)));
display(strcat('pII:', num2str(p_II)));

display(strcat('gEE:', num2str(g_syn_AMPA_on_E)));
display(strcat('gEI:', num2str(g_syn_AMPA_on_I)));
display(strcat('gIE:', num2str(g_syn_GABA_on_E)));
display(strcat('gII:', num2str(g_syn_GABA_on_I)));

end


function [R_t, t] = cal_pop_activity(SpkProfile, N_oscs, tSimBegin, dt, tSimEnd)
[rows,cols]=find(SpkProfile==-1);

% Determine the number of spikes %
Spikes = NaN(N_oscs, 1);
Spikes(1,1) = rows(1) - 1;
Spikes(2:N_oscs,1) = rows(2:N_oscs) - rows(1:(N_oscs - 1)) - 1;

[C, I] = max(Spikes);

Spike_times = ones(N_oscs, Spikes(I, 1)).*(-1);

for i = 1:N_oscs
    if (i == 1)
        if (0 < Spikes(1,1))
            Spike_times(1, 1:Spikes(1,1)) = SpkProfile(1:(rows(1)-1), 1)*1e-3;
        end
    else
        if (0 < Spikes(i,1))
            Spike_times(i, 1:Spikes(i,1)) = SpkProfile((rows(i - 1) + 1):(rows(i) - 1), 1)*1e-3;
        end
    end
end

N_col = size([tSimBegin:dt:tSimEnd], 2);
isNeuronActivation = zeros(1, N_col - 2);
R_t = zeros(1, N_col - 2);

for i = 1:N_oscs
    tmp = hist(Spike_times(i, :), [tSimBegin:dt:tSimEnd]*1e-3);
    isNeuronActivation(1, :) = tmp(1, 2:(N_col - 1));
    isNeuronActivation(1, (0 < isNeuronActivation(1, :))) = 1;
    
    R_t = R_t + isNeuronActivation;
end

R_t = R_t/(N_oscs*dt*1e-3);

tmp = [tSimBegin:dt:tSimEnd];
t = tmp(1, 2:(N_col - 1));

end

function val = isCorrectCheckBlobs(n_E_spk, n_I_spk, i_t_E_blobs_begins, i_t_I_blobs_begins, i_t_E_blobs_ends, i_t_I_blobs_ends)
val = 1;

% Check E
n_e_spks = size(n_E_spk, 2);
n_e = size(i_t_E_blobs_begins, 2);

% Check E before the beginning
i_begin = i_t_E_blobs_begins(1, 1);
for i = 1:(i_begin - 1)
    if (n_E_spk(1,i) ~= 0)
        display('Fail 1');
        val = 0;
        return;
    end
end

% Check E in between
for i = 1:(n_e - 1)
    i_begin = i_t_E_blobs_ends(1, i);
    i_end = i_t_E_blobs_begins(1, i + 1);
    for j = (i_begin + 1):(i_end - 1)
        if (n_E_spk(1,j) ~= 0)
            display('Fail 2');
            val = 0;
            return;
        end
    end
end

% Check E after the end
i_end = i_t_E_blobs_ends(1, n_e);
for i = (i_end + 1):n_e_spks
    if (n_E_spk(1,i) ~= 0)
        display('Fail 3');
        val = 0;
        return;
    end
end

% Check I
n_i_spks = size(n_I_spk, 2);
n_i = size(i_t_I_blobs_begins, 2);

% Check I before the beginning
i_begin = i_t_I_blobs_begins(1, 1);
for i = 1:(i_begin - 1)
    if (n_I_spk(1,i) ~= 0)
        display('Fail 4');
        val = 0;
        return;
    end
end

% Check I in between
for i = 1:(n_i - 1)
    i_begin = i_t_I_blobs_ends(1, i);
    i_end = i_t_I_blobs_begins(1, i + 1);
    for j = (i_begin + 1):(i_end - 1)
        if (n_I_spk(1,j) ~= 0)
            display('Fail 5');
            val = 0;
            return;
        end
    end
end

% Check I after the end
i_end = i_t_I_blobs_ends(1, n_i);
for i = (i_end + 1):n_i_spks
    if (n_I_spk(1,i) ~= 0)
        display('Fail 6');
        val = 0;
        return;
    end
end
end

