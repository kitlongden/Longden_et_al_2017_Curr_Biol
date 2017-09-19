%% 

calc_flag = 1;

N_cells = 10;
N_trials = 10;
N_tf = 9;
N_samples = 20000;
N_ms = 2000;
grey_duration = 0.5; % 0.25 before and 0.25 after
stim_duration = 0.5;
samp_freq = 20000;
    
% Temporal frequencies are: 
tf = [0 1 2 4 7 10 13 17 25];

% Data format is data(Outer N_tf, Inner N_tf, N_samples)
datdir{1} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-1/bindata/''';
datdir{2} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-2/bindata/''';
datdir{3} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-4/bindata/''';
datdir{4} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-6/bindata/''';
datdir{5} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-8/bindata/''';
datdir{6} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-10/bindata/''';
datdir{7} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-1/bindata/''';
datdir{8} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-5/bindata/''';
datdir{9} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-7/bindata/''';
datdir{10}= '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-9/bindata/''';

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
    
    % Data format is (N_cells, Outer N_tf, Inner N_tf, N_samples)
    r_burst_1 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_burst_2 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_burst_3 = zeros(N_cells, N_tf, N_tf, N_trials);
    
    r_trans_1 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_trans_2 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_trans_3 = zeros(N_cells, N_tf, N_tf, N_trials);
    
    r_stead_1 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_stead_2 = zeros(N_cells, N_tf, N_tf, N_trials);
    r_stead_3 = zeros(N_cells, N_tf, N_tf, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    gw = gausswin(200,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_tf,N_tf,N_samples);
    afbarray_2 = zeros(N_tf,N_tf,N_samples);
    afbarray_3 = zeros(N_tf,N_tf,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_tf, N_tf, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(temporal frequency, temporal frequency, number of samples).
            for ntf1 = 1:N_tf
                for ntf2 = 1:N_tf
                    id{ntf1,ntf2,nt} = (find(data(ntf1,ntf2,:) > 0))';
                    id{ntf1,ntf2,nt} = id{ntf1,ntf2,nt} / 20; % conversion to ms
                end
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_motion_parallax(id, N_tf, N_tf, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for ntf1 = 1:N_tf
            for ntf2 = 1:N_tf
                for nt = 1:N_trials
                    % Total response
                    r_burst_1(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,samp_1:samp_2) == 1))/samp_duration;
                    r_burst_2(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,samp_1:samp_2) == 2))/samp_duration;
                    r_burst_3(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,samp_1:samp_2) >= 3))/samp_duration;
                    % Transient
                    r_trans_1(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,trans_1:trans_2) == 1))/trans_duration;
                    r_trans_2(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,trans_1:trans_2) == 2))/trans_duration;
                    r_trans_3(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,trans_1:trans_2) >= 3))/trans_duration;
                    % Steady state
                    r_stead_1(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,stead_1:stead_2) == 1))/stead_duration;
                    r_stead_2(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,stead_1:stead_2) == 2))/stead_duration;
                    r_stead_3(nc,ntf1,ntf2,nt) = length(find(burst_array(ntf1,ntf2,nt,stead_1:stead_2) >= 3))/stead_duration;
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
        for ntf1 = 1:N_tf
            for ntf2 = 1:N_tf
                fbarray_1{nc}(ntf1,ntf2,:) = filtfilt(gw,1,barray_1{nc}(ntf1,ntf2,:));
                fbarray_2{nc}(ntf1,ntf2,:) = filtfilt(gw,1,barray_2{nc}(ntf1,ntf2,:));
                fbarray_3{nc}(ntf1,ntf2,:) = filtfilt(gw,1,barray_3{nc}(ntf1,ntf2,:));
                afbarray_1(ntf1,ntf2,:) = afbarray_1(ntf1,ntf2,:) + fbarray_1{nc}(ntf1,ntf2,:);
                afbarray_2(ntf1,ntf2,:) = afbarray_2(ntf1,ntf2,:) + fbarray_2{nc}(ntf1,ntf2,:);
                afbarray_3(ntf1,ntf2,:) = afbarray_3(ntf1,ntf2,:) + fbarray_3{nc}(ntf1,ntf2,:);
            end
        end
    end
    
    afbarray_1 = afbarray_1 / N_cells;
    afbarray_2 = afbarray_2 / N_cells;
    afbarray_3 = afbarray_3 / N_cells;
    
    disp('    save r_burst data...')
    save r_burst_data_VT1 r_burst* burst_counts barray_* fbarray_* afbarray_* r_trans* r_stead*
end

% PSTHs of mean responses.
% First row of panels is single spikes
% Second row of panels is 2 spike bursts
% Third row of panels is 3+ spike bursts
t = (1:N_samples)/samp_freq;
tf2col = zeros(N_tf,3);
tf2col(:,1) = (1:N_tf)/N_tf;
tf2col(:,3) = (N_tf:-1:1)/N_tf;
lgc = 0.5*[1 1 1];
lgc2 = 0.75*[1 1 1];
tf1col = zeros(N_tf,3);
tf1col(:,1) = (1:N_tf)/N_tf;
tf1col(:,3) = (N_tf:-1:1)/N_tf;

load r_burst_data_VT1.mat afbarray_*
h = figure(5);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
nx = N_tf;
ylim = 1100;
dy = 120;
dyht = num2str(dy/2);

for ntf1 = 1:N_tf
    
    for bn = 1:3
        subplot(3,nx,ntf1 + (bn-1)*nx)
        hold off
        patch([0.25 0.3 0.3 0.25], [0 0 ylim ylim], lgc2, 'edgecolor', lgc2)
        hold on
        for ntf2 = 1:N_tf
            eval(['data = squeeze(afbarray_' num2str(bn) '(ntf1,ntf2,:));'])
            plot(t, data + (ntf2-1)*dy, 'color', tf2col(ntf2,:))
            text(1.01,(ntf2-1)*dy, [num2str(tf(ntf2)) ' Hz'], 'color', tf2col(ntf2,:))
        end
        line([0.25 0.25], [0 ylim], 'color', lgc)
        line([0.75 0.75], [0 ylim], 'color', lgc)
        axis([0 1 0 ylim])
        box off
        xlabel('Time (s)')
        set(gca,'xtick',[0 0.25 0.75 1])
        if (mod(ntf1 + (bn-1)*nx, nx) == 1)
            set(gca,'ytick',0:dy/2:ylim-dy/2, 'yticklabel',{'0', dyht, '0', dyht, '0', dyht, '0', dyht, '0', dyht, '0', dyht, '0', dyht, '0', dyht, '0', dyht} )
            ylabel('Rate (Hz)')
        else
            set(gca,'ytick',0:dy/2:ylim-dy/2, 'yticklabel',{'', '', '', '','', '','', '','', '','', '','', '','', '','', ''} )
        end
        if (tf(ntf1) < 10)
            text(0.05, ylim+50, ['Centre'; num2str(tf(ntf1)) ' Hz  '])
        else
            text(0.05, ylim+50, ['Centre'; num2str(tf(ntf1)) ' Hz '])
        end
        text(0.88, ylim-80, 'Surround', 'color', 'r')
        if (bn < 3)
            text(0.78, ylim+50, [num2str(bn) ' spb'], 'fontweight', 'bold')
        else
            text(0.78, ylim+50, [num2str(bn) '+ spb'], 'fontweight', 'bold')
        end
    end
end

% TF x TF TUNING FOR TRANSIENT RESPONSES
load r_burst_data_VT1.mat  r_trans* r_stead*

h = figure(6);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
nx = 4;
nx2 = 8;

for i = 1:3
    id = num2str(i);
    eval(['mtrans_' id ' = mean(r_trans_' id ' ,4);'])
    eval(['mmtrans_' id ' = squeeze(mean(mtrans_' id '));'])
    eval(['mstead_' id ' = mean(r_stead_' id ' ,4);'])
    eval(['mmstead_' id ' = squeeze(mean(mstead_' id '));'])
end

for i = 1:3
    subplot(3,nx,(i-1)*nx+1)
    eval(['surf(tf,tf,mmtrans_' num2str(i) ')'])
    axis([0 25 0 25 0 150])
    view(2)
    colormap bone
    colorbar
    if (i<3)
        eval(['title(''Transient: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Transient: ' num2str(i) '+ spb'' )'])
    end
    xlabel('Centre temporal frequency (Hz)')
    ylabel('Surround temporal frequency (Hz)')
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+3)
    hold off
    data2 = zeros(1,N_tf);
    for ntf1 = 1:4 %N_tf
        eval(['data = mmtrans_' num2str(i) '(' num2str(ntf1) ',:);'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf1,:), 'markerfacecolor', tf2col(ntf1,:))
        hold on
        text(26,ntf1*8,[num2str(tf(ntf1)) ' Hz'], 'color', tf2col(ntf1,:))
    end
    data2 = data2/4;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Centre temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Transient: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Transient: ' num2str(i) '+ spb'' )'])
    end
    text(26,40,'Surround', 'color', 'b')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+4)
    hold off
    data2 = zeros(1,N_tf);
    for ntf1 = 5:N_tf
        eval(['data = mmtrans_' num2str(i) '(' num2str(ntf1) ',:);'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf1,:), 'markerfacecolor', tf2col(ntf1,:))
        hold on
        text(26,ntf1*8,[num2str(tf(ntf1)) ' Hz'], 'color', tf2col(ntf1,:))
    end
    data2 = data2/5;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Centre temporal frequency (Hz)')
    % ylabel('Rate (Hz)')
    set(gca,'ytick',0:10:90,'yticklabel',{'','','','','','','','','','',})
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Transient: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Transient: ' num2str(i) '+ spb'' )'])
    end
    text(26,80,'Surround', 'color', 'r')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+6)
    hold off
    data2 = zeros(N_tf,1);
    for ntf2 = 1:4 %N_tf
        eval(['data = mmtrans_' num2str(i) '(:,' num2str(ntf2) ');'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf2,:), 'markerfacecolor', tf2col(ntf2,:))
        hold on
        text(26,ntf2*8,[num2str(tf(ntf2)) ' Hz'], 'color', tf2col(ntf2,:))
    end
    data2 = data2/4;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Surround temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Transient: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Transient: ' num2str(i) '+ spb'' )'])
    end
    text(26,40,'Centre', 'color', 'b')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+7)
    hold off
    data2 = zeros(N_tf,1);
    for ntf2 = 5:N_tf
        eval(['data = mmtrans_' num2str(i) '(:,' num2str(ntf2) ');'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf2,:), 'markerfacecolor', tf2col(ntf2,:))
        hold on
        text(26,ntf2*8,[num2str(tf(ntf2)) ' Hz'], 'color', tf2col(ntf2,:))
    end
    data2 = data2/5;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Surround temporal frequency (Hz)')
    % ylabel('Rate (Hz)')
    set(gca,'ytick',0:10:90,'yticklabel',{'','','','','','','','','','',})
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Transient: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Transient: ' num2str(i) '+ spb'' )'])
    end
    text(26,80,'Centre', 'color', 'r')
    box off
end

% TF x TF TUNING FOR STEADY STATE RESPONSES

h = figure(7);
set(h,'Position', [60, 20, 2000, 1000])
set(h,'color',[1 1 1])
nx = 4;
nx2 = 8

for i = 1:3
    subplot(3,nx,(i-1)*nx+1)
    eval(['surf(tf,tf,mmstead_' num2str(i) ')'])
    axis([0 25 0 25 0 150])
    view(2)
    colormap bone
    colorbar
    if (i<3)
        eval(['title(''Steady state: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Steady state: ' num2str(i) '+ spb'' )'])
    end
    xlabel('Centre temporal frequency (Hz)')
    ylabel('Surround temporal frequency (Hz)')
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+3)
    hold off
    data2 = zeros(1,N_tf);
    for ntf1 = 1:4 %N_tf
        eval(['data = mmstead_' num2str(i) '(' num2str(ntf1) ',:);'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf1,:), 'markerfacecolor', tf2col(ntf1,:))
        hold on
        text(26,ntf1*8,[num2str(tf(ntf1)) ' Hz'], 'color', tf2col(ntf1,:))
    end
    data2 = data2/4;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Centre temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Steady state: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Steady state: ' num2str(i) '+ spb'' )'])
    end
    text(26,40,'Surround', 'color', 'b')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+4)
    hold off
    data2 = zeros(1,N_tf);
    for ntf1 = 5:N_tf
        eval(['data = mmstead_' num2str(i) '(' num2str(ntf1) ',:);'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf1,:), 'markerfacecolor', tf2col(ntf1,:))
        hold on
        text(26,ntf1*8,[num2str(tf(ntf1)) ' Hz'], 'color', tf2col(ntf1,:))
    end
    data2 = data2/5;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Centre temporal frequency (Hz)')
    % ylabel('Rate (Hz)')
    set(gca,'ytick',0:10:90,'yticklabel',{'','','','','','','','','','',})
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Steady state: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Steady state: ' num2str(i) '+ spb'' )'])
    end
    text(26,80,'Surround', 'color', 'r')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+6)
    hold off
    data2 = zeros(N_tf,1);
    for ntf2 = 1:4 %N_tf
        eval(['data = mmstead_' num2str(i) '(:,' num2str(ntf2) ');'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf2,:), 'markerfacecolor', tf2col(ntf2,:))
        hold on
        text(26,ntf2*8,[num2str(tf(ntf2)) ' Hz'], 'color', tf2col(ntf2,:))
    end
    data2 = data2/4;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Surround temporal frequency (Hz)')
    ylabel('Rate (Hz)')
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Steady state: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Steady state: ' num2str(i) '+ spb'' )'])
    end
    text(26,40,'Centre', 'color', 'b')
    box off
end

for i = 1:3
    subplot(3,nx2,(i-1)*nx2+7)
    hold off
    data2 = zeros(N_tf,1);
    for ntf2 = 5:N_tf
        eval(['data = mmstead_' num2str(i) '(:,' num2str(ntf2) ');'])
        data2 = data2 + data;
        plot(tf, data, 'o-', 'color', tf2col(ntf2,:), 'markerfacecolor', tf2col(ntf2,:))
        hold on
        text(26,ntf2*8,[num2str(tf(ntf2)) ' Hz'], 'color', tf2col(ntf2,:))
    end
    data2 = data2/5;
    plot(tf, data2, 'ko-', 'markerfacecolor', 'k', 'linewidth', 2)
    axis([0 26 0 90])
    xlabel('Surround temporal frequency (Hz)')
    % ylabel('Rate (Hz)')
    set(gca,'ytick',0:10:90,'yticklabel',{'','','','','','','','','','',})
    set(gca,'xtick', tf)
    if (i<3)
        eval(['title(''Steady state: ' num2str(i) ' spb'' )'])
    else
        eval(['title(''Steady state: ' num2str(i) '+ spb'' )'])
    end
    text(26,80,'Centre', 'color', 'r')
    box off
end























