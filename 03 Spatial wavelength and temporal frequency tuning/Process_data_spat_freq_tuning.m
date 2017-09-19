%% 

calc_flag = 1;

N_cells = 10;
N_trials = 25;
N_tf = 8;
N_sw = 6;
N_samples = 20000;
N_ms = 2000;
grey_duration = 0.5; % 0.25 before and 0.25 after
stim_duration = 0.5;
samp_freq = 20000;
    
tf = [0.1 0.5 2 4 7 10 17 25];
sw = [5.6160   11.1262   13.8113   16.4361   23.8694   30.5406];

% Temporal frequencies are: 
% 0.1 0.5 2 4 7 10 17 25
% 
% Spatial wavelengths at the top of the screen:  
% 5.6160   11.1262   13.8113   16.4361   23.8694   30.5406

datdir{1} =   '''../Data/Spatial wavelength and temporal frequency tuning/06 May 2014-3/bindata/''';
datdir{2} =   '''../Data/Spatial wavelength and temporal frequency tuning/14 May 2014-7/bindata/''';
datdir{3} =   '''../Data/Spatial wavelength and temporal frequency tuning/15 May 2014-2/bindata/''';
datdir{4} =   '''../Data/Spatial wavelength and temporal frequency tuning/21 Apr 2014-3/bindata/''';
datdir{5} =   '''../Data/Spatial wavelength and temporal frequency tuning/22 Apr 2014-3/bindata/''';
datdir{6} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-10/bindata/''';
datdir{7} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-2/bindata/''';
datdir{8} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-7/bindata/''';
datdir{9} =   '''../Data/Spatial wavelength and temporal frequency tuning/24 Apr 2014-2/bindata/''';
datdir{10} =  '''../Data/Spatial wavelength and temporal frequency tuning/24 Apr 2014-5/bindata/''';

if (calc_flag == 1)
    
    samp_1 = 5001;
    samp_2 = 15000;
    samp_duration = (samp_2-samp_1)/samp_freq;
    trans_1 = 5001;
    trans_2 = 6000;
    trans_duration = (trans_2-trans_1)/samp_freq;
    stead_1 = 6001;
    stead_2 = 15000;
    stead_duration = (stead_2-stead_1)/samp_freq;
    spont_1 = 1;
    spont_2 = 5000;
    spont_3 = 15001;
    spont_4 = 20000;
    spont_duration = ((spont_2-spont_1)+(spont_4-spont_3))/samp_freq;
    
    r_burst_1 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_burst_2 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_burst_3 = zeros(N_cells, N_sw, N_tf, N_trials);
    
    r_trans_1 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_trans_2 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_trans_3 = zeros(N_cells, N_sw, N_tf, N_trials);
    
    r_stead_1 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_stead_2 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_stead_3 = zeros(N_cells, N_sw, N_tf, N_trials);
    
    r_spont_1 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_spont_2 = zeros(N_cells, N_sw, N_tf, N_trials);
    r_spont_3 = zeros(N_cells, N_sw, N_tf, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    gw = gausswin(200,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_sw,N_tf,N_samples);
    afbarray_2 = zeros(N_sw,N_tf,N_samples);
    afbarray_3 = zeros(N_sw,N_tf,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_sw, N_tf, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(spatial frequency, temporal frequency, number of samples).
            for nsw = 1:N_sw
                for ntf = 1:N_tf
                    id{nsw,ntf,nt} = (find(data(nsw,ntf,:) > 0))';
                    id{nsw,ntf,nt} = id{nsw,ntf,nt} / 20; % conversion to ms
                end
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_spat_freq(id, N_sw, N_tf, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for nsw = 1:N_sw
            for ntf = 1:N_tf
                for nt = 1:N_trials
                    % Total response
                    r_burst_1(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,samp_1:samp_2) == 1))/samp_duration;
                    r_burst_2(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,samp_1:samp_2) == 2))/samp_duration;
                    r_burst_3(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,samp_1:samp_2) >= 3))/samp_duration;
                    % Transient
                    r_trans_1(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,trans_1:trans_2) == 1))/trans_duration;
                    r_trans_2(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,trans_1:trans_2) == 2))/trans_duration;
                    r_trans_3(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,trans_1:trans_2) >= 3))/trans_duration;
                    % Steady state
                    r_stead_1(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,stead_1:stead_2) == 1))/stead_duration;
                    r_stead_2(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,stead_1:stead_2) == 2))/stead_duration;
                    r_stead_3(nc,nsw,ntf,nt) = length(find(burst_array(nsw,ntf,nt,stead_1:stead_2) >= 3))/stead_duration;
                    % Spontaneous
                    r_spont_1(nc,nsw,ntf,nt) = (length(find(burst_array(nsw,ntf,nt,spont_1:spont_2) == 1))+length(find(burst_array(nsw,ntf,nt,spont_3:spont_4) == 1)))/spont_duration;
                    r_spont_2(nc,nsw,ntf,nt) = (length(find(burst_array(nsw,ntf,nt,spont_1:spont_2) == 2))+length(find(burst_array(nsw,ntf,nt,spont_3:spont_4) == 2)))/spont_duration;
                    r_spont_3(nc,nsw,ntf,nt) = (length(find(burst_array(nsw,ntf,nt,spont_1:spont_2) >= 3))+length(find(burst_array(nsw,ntf,nt,spont_3:spont_4) >= 3)))/spont_duration;
                end
            end
        end
        
        disp('   Calculate ISI distributions...')
        for nsb = 1:10
            burst_counts(nc,nsb) = length(find(burst_array == nsb));
        end
        
        disp('   Calculate mean response traces...')
        burst_array_1 = (burst_array == 1)*samp_freq;
        burst_array_2 = (burst_array == 2)*samp_freq;
        burst_array_3 = (burst_array >= 3)*samp_freq;
        barray_1{nc} = squeeze(mean(burst_array_1,3));
        barray_2{nc} = squeeze(mean(burst_array_2,3));
        barray_3{nc} = squeeze(mean(burst_array_3,3));
        
        disp('   Filter mean response traces...')
        for nsw = 1:N_sw
            for ntf = 1:N_tf
                fbarray_1{nc}(nsw,ntf,:) = filtfilt(gw,1,barray_1{nc}(nsw,ntf,:));
                fbarray_2{nc}(nsw,ntf,:) = filtfilt(gw,1,barray_2{nc}(nsw,ntf,:));
                fbarray_3{nc}(nsw,ntf,:) = filtfilt(gw,1,barray_3{nc}(nsw,ntf,:));
                afbarray_1(nsw,ntf,:) = afbarray_1(nsw,ntf,:) + fbarray_1{nc}(nsw,ntf,:);
                afbarray_2(nsw,ntf,:) = afbarray_2(nsw,ntf,:) + fbarray_2{nc}(nsw,ntf,:);
                afbarray_3(nsw,ntf,:) = afbarray_3(nsw,ntf,:) + fbarray_3{nc}(nsw,ntf,:);
            end
        end
    end
    
    afbarray_1 = afbarray_1 / N_cells;
    afbarray_2 = afbarray_2 / N_cells;
    afbarray_3 = afbarray_3 / N_cells;
    
    disp('    save r_burst data...')
    save r_burst_data r_burst* burst_counts barray_* fbarray_* afbarray_* r_trans* r_stead* r_spont*
end


% MEAN RESPONSES FOR ALL SPATIAL & TEMPORAL FREQUENCIES & EVERY CELL
t = (1:N_samples)/samp_freq;
tfcol = zeros(N_tf,3);
tfcol(:,1) = (1:N_tf)/N_tf;
tfcol(:,3) = (N_tf:-1:1)/N_tf;
lgc = 0.5*[1 1 1];
lgc2 = 0.75*[1 1 1];
swcol = zeros(N_sw,3);
swcol(:,1) = (1:N_sw)/N_sw;
swcol(:,3) = (N_sw:-1:1)/N_sw;

load r_burst_data.mat afbarray_*
h = figure(N_cells+1+1);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
for nsw = 2:N_sw
    
    for bn = 1:3
        subplot(3,5,(nsw-1) + (bn-1)*5)
        hold off
        patch([0.25 0.3 0.3 0.25], [0 0 1000 1000], lgc2, 'edgecolor', lgc2)
        hold on
        for ntf = 1:N_tf
            eval(['data = squeeze(afbarray_' num2str(bn) '(nsw,ntf,:));'])
            plot(t, data + (ntf-1)*50, 'color', tfcol(ntf,:))
            text(1.01,(ntf-1)*50, [num2str(tf(ntf)) ' Hz'], 'color', tfcol(ntf,:))
        end
        line([0.25 0.25], [0 1000], 'color', lgc)
        line([0.75 0.75], [0 1000], 'color', lgc)
        axis([0 1 0 425])
        box off
        xlabel('Time (s)')
        ylabel('Rate (Hz)')
        set(gca,'xtick',[0 0.25 0.75 1])
        set(gca,'ytick',0:25:375, 'yticklabel',{'0', '25', '0', '25', '0', '25', '0', '25', '0', '25', '0', '25', '0', '25', '0', '25'} )
        text(0.05, 400, [num2str(round(sw(nsw)*10)/10) ' ^o'])
        if (bn < 3)
            text(0.85, 400, [num2str(bn) ' spb'], 'fontweight', 'bold')
        else
            text(0.85, 400, [num2str(bn) '+ spb'], 'fontweight', 'bold')
        end
    end
end

h = figure(N_cells+1+1+1);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
for ntf = 1:N_tf    
    for bn = 1:3
        subplot(3,N_tf,ntf + (bn-1)*N_tf)
        hold off
        patch([0.25 0.3 0.3 0.25], [0 0 1000 1000], lgc2, 'edgecolor', lgc2)
        hold on
        for nsw = 1:N_sw
            eval(['data = squeeze(afbarray_' num2str(bn) '(nsw,ntf,:));'])
            plot(t, data + (nsw)*50, 'color', swcol(nsw,:))
            text(1.01,(nsw)*50, [num2str(round(sw(nsw))) '^o'], 'color', swcol(nsw,:))
        end
        line([0.25 0.25], [0 1000], 'color', lgc)
        line([0.75 0.75], [0 1000], 'color', lgc)
        axis([0 1 0 375])
        box off
        xlabel('Time (s)')
        ylabel('Rate (Hz)')
        set(gca,'xtick',[0 0.25 0.75 1])
        set(gca,'ytick',50:25:325, 'yticklabel',{'0', '25', '0', '25', '0', '25', '0', '25', '0', '25', '0', '25'} )
        text(0.05, 350, [num2str(tf(ntf)) ' Hz'])
        if (bn < 3)
            text(0.8, 350, [num2str(bn) ' spb'], 'fontweight', 'bold')
        else
            text(0.8, 350, [num2str(bn) '+ spb'], 'fontweight', 'bold')
        end
    end
end

% TF & SW TUNING FOR TRANSIENT AND STEADY STATE RESPONSES
h = figure(N_cells+1+1+1+1);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
nx = 8;
load r_burst_data.mat r_trans* r_stead* 

for i = 1:3
    id = num2str(i);
    eval(['mtrans_' id ' = mean(r_trans_' id ' ,4);'])
    eval(['mmtrans_' id ' = squeeze(mean(mtrans_' id '));'])
    eval(['mstead_' id ' = mean(r_stead_' id ' ,4);'])
    eval(['mmstead_' id ' = squeeze(mean(mstead_' id '));'])
end

for i = 1:3
    subplot(3,nx,(i-1)*nx+1)
    eval(['surf(tf,sw,mmtrans_' num2str(i) ')'])
    axis([0 25 5 31 0 50])
    view(2)
    colormap bone
    colorbar
    eval(['title(''Transient: ' num2str(i) ' spikes per burst'' )']) 
    xlabel('Temporal frequency (Hz)')
    ylabel('Spatial frequency (^o)')
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+2)
    hold off
    eval(['data6 = mmtrans_' num2str(i) '(6,:);'])
    plot(tf, data6, 'o-', 'color', swcol(6,:), 'markerfacecolor', swcol(6,:))
    hold on
    eval(['data5 = mmtrans_' num2str(i) '(5,:);'])
    plot(tf, data5, 'o-', 'color', swcol(5,:), 'markerfacecolor', swcol(5,:))
    eval(['data4 = mmtrans_' num2str(i) '(4,:);'])
    plot(tf, data4, 'o-', 'color', swcol(4,:), 'markerfacecolor', swcol(4,:))
    eval(['data3 = mmtrans_' num2str(i) '(3,:);'])
    plot(tf, data3, 'o-', 'color', swcol(3,:), 'markerfacecolor', swcol(3,:))
    eval(['data2 = mmtrans_' num2str(i) '(2,:);'])
    plot(tf, data2, 'o-', 'color', swcol(2,:), 'markerfacecolor', swcol(2,:))
    eval(['data1 = mmtrans_' num2str(i) '(1,:);'])
    plot(tf, data1, 'o-', 'color', swcol(1,:), 'markerfacecolor', swcol(1,:))
    data456 = (data4+data5+data6)/3;
    plot(tf, data456, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 26 0 20])
    xlabel('Temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    eval(['title(''Transient: ' num2str(i) ' spikes per burst'' )'])
    text(26,16,' Mean')
    text(26,14,[num2str(round(sw(6))) ' ^o'], 'color', swcol(6,:))
    text(26,12,[num2str(round(sw(5))) ' ^o'], 'color', swcol(5,:))
    text(26,10,[num2str(round(sw(4))) ' ^o'], 'color', swcol(4,:))
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+3)
    hold off
    eval(['data8 = mmtrans_' num2str(i) '(:,8);'])
    plot(sw, data8, 'o-', 'color', tfcol(8,:), 'markerfacecolor', tfcol(8,:))
    hold on
    eval(['data7 = mmtrans_' num2str(i) '(:,7);'])
    plot(sw, data7, 'o-', 'color', tfcol(7,:), 'markerfacecolor', tfcol(7,:))
    eval(['data6 = mmtrans_' num2str(i) '(:,6);'])
    plot(sw, data6, 'o-', 'color', tfcol(6,:), 'markerfacecolor', tfcol(6,:))
    eval(['data5 = mmtrans_' num2str(i) '(:,5);'])
    plot(sw, data5, 'o-', 'color', tfcol(5,:), 'markerfacecolor', tfcol(5,:))
    data5678 = (data5+data6+data7+data8)/4;
    plot(sw, data5678, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 31 0 20])
    xlabel('Spatial frequency (^o)')
    eval(['title(''Transient: ' num2str(i) ' spikes per burst'' )'])
    text(31,16,' Mean')
    text(31,14,[num2str(tf(8)) ' Hz'], 'color', tfcol(8,:))
    text(31,12,[num2str(tf(7)) ' Hz'], 'color', tfcol(7,:))
    text(31,10,[num2str(tf(6)) ' Hz'], 'color', tfcol(6,:))
    text(31,8,[num2str(tf(5)) ' Hz'], 'color', tfcol(5,:))
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+4)
    hold off
    eval(['data4 = mmtrans_' num2str(i) '(:,4);'])
    plot(sw, data4, 'o-', 'color', tfcol(4,:), 'markerfacecolor', tfcol(4,:))
    hold on
    eval(['data3 = mmtrans_' num2str(i) '(:,3);'])
    plot(sw, data3, 'o-', 'color', tfcol(3,:), 'markerfacecolor', tfcol(3,:))
    eval(['data2 = mmtrans_' num2str(i) '(:,2);'])
    plot(sw, data2, 'o-', 'color', tfcol(2,:), 'markerfacecolor', tfcol(2,:))
    eval(['data1 = mmtrans_' num2str(i) '(:,1);'])
    plot(sw, data1, 'o-', 'color', tfcol(1,:), 'markerfacecolor', tfcol(1,:))
    data1234 = (data4+data3+data2+data1)/4;
    plot(sw, data1234, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 31 0 20])
    xlabel('Spatial frequency (^o)')
    eval(['title(''Transient: ' num2str(i) ' spikes per burst'' )'])
    text(31,16,' Mean')
    text(31,14,[num2str(tf(4)) ' Hz'], 'color', tfcol(4,:))
    text(31,12,[num2str(tf(3)) ' Hz'], 'color', tfcol(3,:))
    text(31,10,[num2str(tf(2)) ' Hz'], 'color', tfcol(2,:))
    text(31,8,[num2str(tf(1)) ' Hz'], 'color', tfcol(1,:))
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+5)
    eval(['surf(tf,sw,mmstead_' num2str(i) ')'])
    axis([0 25 5 31 0 50])
    view(2)
    colormap bone
    colorbar
    eval(['title(''Steady state: ' num2str(i) ' spikes per burst'' )']) 
    xlabel('Temporal frequency (Hz)')
    ylabel('Spatial frequency (^o)')
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+6)
    hold off
    eval(['data6 = mmstead_' num2str(i) '(6,:);'])
    plot(tf, data6, 'o-', 'color', swcol(6,:), 'markerfacecolor', swcol(6,:))
    hold on
    eval(['data5 = mmstead_' num2str(i) '(5,:);'])
    plot(tf, data5, 'o-', 'color', swcol(5,:), 'markerfacecolor', swcol(5,:))
    eval(['data4 = mmstead_' num2str(i) '(4,:);'])
    plot(tf, data4, 'o-', 'color', swcol(4,:), 'markerfacecolor', swcol(4,:))
    data456 = (data4+data5+data6)/3;
    plot(tf, data456, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 26 0 30])
    xlabel('Temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    eval(['title(''Steady state: ' num2str(i) ' spikes per burst'' )'])
    text(26,16,' Mean')
    text(26,14,[num2str(round(sw(6))) ' ^o'], 'color', swcol(6,:))
    text(26,12,[num2str(round(sw(5))) ' ^o'], 'color', swcol(5,:))
    text(26,10,[num2str(round(sw(4))) ' ^o'], 'color', swcol(4,:))
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+7)
    hold off
    eval(['data8 = mmstead_' num2str(i) '(:,8);'])
    plot(sw, data8, 'o-', 'color', tfcol(8,:), 'markerfacecolor', tfcol(8,:))
    hold on
    eval(['data7 = mmstead_' num2str(i) '(:,7);'])
    plot(sw, data7, 'o-', 'color', tfcol(7,:), 'markerfacecolor', tfcol(7,:))
    eval(['data6 = mmstead_' num2str(i) '(:,6);'])
    plot(sw, data6, 'o-', 'color', tfcol(6,:), 'markerfacecolor', tfcol(6,:))
    eval(['data5 = mmstead_' num2str(i) '(:,5);'])
    plot(sw, data5, 'o-', 'color', tfcol(5,:), 'markerfacecolor', tfcol(5,:))
    data5678 = (data5+data6+data7+data8)/4;
    plot(sw, data5678, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 31 0 30])
    xlabel('Spatial frequency (^o)')
    eval(['title(''Steady state: ' num2str(i) ' spikes per burst'' )'])
    text(31,16,' Mean')
    text(31,14,[num2str(tf(8)) ' Hz'], 'color', tfcol(8,:))
    text(31,12,[num2str(tf(7)) ' Hz'], 'color', tfcol(7,:))
    text(31,10,[num2str(tf(6)) ' Hz'], 'color', tfcol(6,:))
    text(31,8,[num2str(tf(5)) ' Hz'], 'color', tfcol(5,:))
end
for i = 1:3
    subplot(3,nx,(i-1)*nx+8)
    hold off
    eval(['data4 = mmstead_' num2str(i) '(:,4);'])
    plot(sw, data4, 'o-', 'color', tfcol(4,:), 'markerfacecolor', tfcol(4,:))
    hold on
    eval(['data3 = mmstead_' num2str(i) '(:,3);'])
    plot(sw, data3, 'o-', 'color', tfcol(3,:), 'markerfacecolor', tfcol(3,:))
    eval(['data2 = mmstead_' num2str(i) '(:,2);'])
    plot(sw, data2, 'o-', 'color', tfcol(2,:), 'markerfacecolor', tfcol(2,:))
    eval(['data1 = mmstead_' num2str(i) '(:,1);'])
    plot(sw, data1, 'o-', 'color', tfcol(1,:), 'markerfacecolor', tfcol(1,:))
    data1234 = (data4+data3+data2+data1)/4;
    plot(sw, data1234, 'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 2)
    box off
    axis([0 31 0 30])
    xlabel('Spatial frequency (^o)')
    eval(['title(''Steady state: ' num2str(i) ' spikes per burst'' )'])
    text(31,16,' Mean')
    text(31,14,[num2str(tf(4)) ' Hz'], 'color', tfcol(4,:))
    text(31,12,[num2str(tf(3)) ' Hz'], 'color', tfcol(3,:))
    text(31,10,[num2str(tf(2)) ' Hz'], 'color', tfcol(2,:))
    text(31,8,[num2str(tf(1)) ' Hz'], 'color', tfcol(1,:))
end


















































