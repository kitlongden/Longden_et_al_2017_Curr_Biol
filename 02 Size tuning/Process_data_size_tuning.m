%% 

% STIMULI
% Pixel radius of circles:  pix_rad=[24,48,72,96,120,144,168,192,216,240];
% Screeheight = 263 mm = 480 pixels; Screendist = 200;
% Top of screen is at elevation 0^o
% Elevation angle of screen = atand(263/200) = 52.7^o.
% At screen centre: elevation=33.3^o, dist=sqrt((263/2)^2 + 200^2)=239
% Circle diameter at screen centre: atand((pix_rad/240)*263/239)
% = [6.2797,12.4120,18.2694,23.7575,28.82,33.4348,37.6068,41.3586,44.723]
% The ring stimuli are made by having a full width circular stimulus, and
% then increasing the diameter of the mask at the centre.
circ_diam = [6.2797,12.4120,18.2694,23.7575,28.82,33.4348,37.6068,41.3586,44.723,47.7372];
circ_diam2 = round(circ_diam*10)/10;
mask_diam = [6.2797,12.4120,18.2694,23.7575,28.82,33.4348,37.6068,41.3586,44.723,47.7372];
mask_diam2 = round(mask_diam*10)/10;

calc_flag = 1;

N_cells = 10;
N_trials = 30;
N_size = 20;
N_diam = 10;
N_samples = 30000;
N_ms = 3000;
grey_duration = 1.0; % 0.5 before and 0.5 after
stim_duration = 0.5;
samp_freq = 20000;

% The data format is as follows:
% data(sizes, samples).

datdir{1} = '''../Data/Size tuning/14 Apr 2014-2/bindata/''';
datdir{2} = '''../Data/Size tuning/17 Mar 2014-2/bindata/''';
datdir{3} = '''../Data/Size tuning/19 Mar 2014-2/bindata/''';
datdir{4} = '''../Data/Size tuning/24 Mar 2014-2/bindata/''';
datdir{5} = '''../Data/Size tuning/24 Mar 2014-3/bindata/''';
datdir{6} = '''../Data/Size tuning/26 Mar 2014-3/bindata/''';
datdir{7} = '''../Data/Size tuning/19 Mar 2014-6/bindata/''';
datdir{8} = '''../Data/Size tuning/21 Mar 2014-2/bindata/''';
datdir{9} = '''../Data/Size tuning/25 Mar 2014-1/bindata/''';
datdir{10}= '''../Data/Size tuning/27 Mar 2014-2/bindata/''';

if (calc_flag == 1)
    
    samp_1 = 10001;
    samp_2 = 20000;
    samp_duration = (samp_2-samp_1)/samp_freq;
    trans_1 = 10001;
    trans_2 = 11000;
    trans_duration = (trans_2-trans_1)/samp_freq;
    stead_1 = 11001;
    stead_2 = 20000;
    stead_duration = (stead_2-stead_1)/samp_freq;
    
    r_burst_1 = zeros(N_cells, N_size, N_trials);
    r_burst_2 = zeros(N_cells, N_size, N_trials);
    r_burst_3 = zeros(N_cells, N_size, N_trials);
    
    r_trans_1 = zeros(N_cells, N_size, N_trials);
    r_trans_2 = zeros(N_cells, N_size, N_trials);
    r_trans_3 = zeros(N_cells, N_size, N_trials);
    
    r_stead_1 = zeros(N_cells, N_size, N_trials);
    r_stead_2 = zeros(N_cells, N_size, N_trials);
    r_stead_3 = zeros(N_cells, N_size, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    % gw = gausswin(200,2.5);
    gw = gausswin(400,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_size,N_samples);
    afbarray_2 = zeros(N_size,N_samples);
    afbarray_3 = zeros(N_size,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_size, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(size, number of samples).
            for ns = 1:N_size
                id{ns,nt} = (find(data(ns,:) > 0))';
                id{ns,nt} = id{ns,nt} / 20; % conversion to ms
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_size_tuning(id, N_size, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for nsz = 1:N_size
            for nt = 1:N_trials
                % Total response
                r_burst_1(nc,nsz,nt) = length(find(burst_array(nsz,nt,samp_1:samp_2) == 1))/samp_duration;
                r_burst_2(nc,nsz,nt) = length(find(burst_array(nsz,nt,samp_1:samp_2) == 2))/samp_duration;
                r_burst_3(nc,nsz,nt) = length(find(burst_array(nsz,nt,samp_1:samp_2) >= 3))/samp_duration;
                % Transient
                r_trans_1(nc,nsz,nt) = length(find(burst_array(nsz,nt,trans_1:trans_2) == 1))/trans_duration;
                r_trans_2(nc,nsz,nt) = length(find(burst_array(nsz,nt,trans_1:trans_2) == 2))/trans_duration;
                r_trans_3(nc,nsz,nt) = length(find(burst_array(nsz,nt,trans_1:trans_2) >= 3))/trans_duration;
                % Steady state
                r_stead_1(nc,nsz,nt) = length(find(burst_array(nsz,nt,stead_1:stead_2) == 1))/stead_duration;
                r_stead_2(nc,nsz,nt) = length(find(burst_array(nsz,nt,stead_1:stead_2) == 2))/stead_duration;
                r_stead_3(nc,nsz,nt) = length(find(burst_array(nsz,nt,stead_1:stead_2) >= 3))/stead_duration;
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
        barray_1{nc} = squeeze(mean(burst_array_1,2));
        barray_2{nc} = squeeze(mean(burst_array_2,2));
        barray_3{nc} = squeeze(mean(burst_array_3,2));
        
        disp('   Filter mean response traces...')
        for nsz = 1:N_size
            fbarray_1{nc}(nsz,:) = filtfilt(gw,1,barray_1{nc}(nsz,:));
            fbarray_2{nc}(nsz,:) = filtfilt(gw,1,barray_2{nc}(nsz,:));
            fbarray_3{nc}(nsz,:) = filtfilt(gw,1,barray_3{nc}(nsz,:));
            afbarray_1(nsz,:) = afbarray_1(nsz,:) + fbarray_1{nc}(nsz,:);
            afbarray_2(nsz,:) = afbarray_2(nsz,:) + fbarray_2{nc}(nsz,:);
            afbarray_3(nsz,:) = afbarray_3(nsz,:) + fbarray_3{nc}(nsz,:);
        end
    end
    
    afbarray_1 = afbarray_1 / N_cells;
    afbarray_2 = afbarray_2 / N_cells;
    afbarray_3 = afbarray_3 / N_cells;
    
    disp('    save r_burst data...')
    save r_burst_data r_burst* burst_counts barray_* fbarray_* afbarray_* r_trans* r_stead*
end

% MEAN RESPONSES FOR ALL SPATIAL & TEMPORAL FREQUENCIES & EVERY CELL
t = (1:N_samples)/samp_freq;
ndcol = zeros(N_diam,3);
ndcol(:,1) = (1:N_diam)/N_diam;
ndcol(:,3) = (N_diam:-1:1)/N_diam;
lgc = 0.5*[1 1 1];
lgc2 = 0.75*[1 1 1];

load r_burst_data.mat afbarray_*
h = figure(N_cells+1+1);
set(h,'Position', [60, 20, 1500, 1000])
set(h,'color',[1 1 1])
nx = 4;
dy = 150;
dyt = num2str(dy/2);
ylim = 1550;

% PSTHs
for bn = 1:3
    subplot(3,nx,1 + (bn-1)*nx)
    hold off
    patch([0.5 0.55 0.55 0.5], [0 0 ylim ylim], lgc2, 'edgecolor', lgc2)
    hold on
    for nd = 1:N_diam
        eval(['data = squeeze(afbarray_' num2str(bn) '(nd,:));'])
        plot(t, data + (nd-1)*dy, 'color', ndcol(nd,:))
        text(1.51,(nd-1)*dy, [num2str(circ_diam2(nd)) ' ^o'], 'color', ndcol(nd,:))
    end
    line([0.5 0.5], [0 ylim], 'color', lgc)
    line([1 1], [0 ylim], 'color', lgc)
    axis([0 1.5 0 ylim])
    box off
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    set(gca,'xtick',[0 0.5 1 1.5])
    set(gca,'ytick',0:dy/2:1450, 'yticklabel',{'0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt} )
    text(0.05, 1500, 'Centre')
    if (bn < 3)
        text(1.025, 1500, [num2str(bn) ' spb'], 'fontweight', 'bold')
    else
        text(1.025, 1500, [num2str(bn) '+ spb'], 'fontweight', 'bold')
    end
    text(1.5,1500,['Centre  ' ; 'diameter'], 'color', 'r')
end
for bn = 1:3
    subplot(3,nx,2 + (bn-1)*nx)
    hold off
    patch([0.5 0.55 0.55 0.5], [0 0 ylim ylim], lgc2, 'edgecolor', lgc2)
    hold on
    for nd = 1:N_diam
        nd_id = N_diam + nd; %N_size - (nd+1);
        eval(['data = squeeze(afbarray_' num2str(bn) '(nd_id,:));'])
        plot(t, data + (nd-1)*dy, 'color', ndcol(nd,:))
        text(1.51,(nd-1)*dy, [num2str(mask_diam2(nd)) ' ^o'], 'color', ndcol(nd,:))
    end
    line([0.5 0.5], [0 ylim], 'color', lgc)
    line([1 1], [0 ylim], 'color', lgc)
    axis([0 1.5 0 ylim])
    box off
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    set(gca,'xtick',[0 0.5 1 1.5])
    set(gca,'ytick',0:dy/2:1450, 'yticklabel',{'0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt} )
    text(0.05, 1500, 'Surround')
    if (bn < 3)
        text(1.025, 1500, [num2str(bn) ' spb'], 'fontweight', 'bold')
    else
        text(1.025, 1500, [num2str(bn) '+ spb'], 'fontweight', 'bold')
    end
    text(1.4,1500,['Centre mask' ; 'diameter   '], 'color', 'r')
end

% TF & SW TUNING FOR TRANSIENT AND STEADY STATE RESPONSES
load r_burst_data.mat r_trans* r_stead* r_burst*
for i = 1:3
    id = num2str(i);
    eval(['mburst_' id ' = mean(r_burst_' id ' ,3);'])
    eval(['mmburst_' id ' = squeeze(mean(mburst_' id '));'])
    eval(['mtrans_' id ' = mean(r_trans_' id ' ,3);'])
    eval(['mmtrans_' id ' = squeeze(mean(mtrans_' id '));'])
    eval(['mstead_' id ' = mean(r_stead_' id ' ,3);'])
    eval(['mmstead_' id ' = squeeze(mean(mstead_' id '));'])
end
mburst_all = mean(r_burst_1+r_burst_2+r_burst_3,3);
mmburst_all = squeeze(mean(mburst_all));
mtrans_all = mean(r_trans_1+r_trans_2+r_trans_3,3);
mmtrans_all = squeeze(mean(mtrans_all));
mstead_all = mean(r_stead_1+r_stead_2+r_stead_3,3);
mmstead_all = squeeze(mean(mstead_all));

bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];
ylim2 = 110;

subplot(3,nx,3)
hold off
for bn = 1:3
    eval(['data = mmtrans_' num2str(bn) ';']);
    plot([0 circ_diam2], [data(20) data(1:10)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmtrans_all;
plot([0 circ_diam2], [data(20) data(1:10)], 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 50 0 ylim2])
xlabel('Centre diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Centre: transient response')

subplot(3,nx,3+nx)
hold off
for bn = 1:3
    eval(['data = mmstead_' num2str(bn) ';']);
    plot([0 circ_diam2], [data(20) data(1:10)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmstead_all;
plot([0 circ_diam2], [data(20) data(1:10)], 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 50 0 ylim2])
xlabel('Centre diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Centre: steady state response')

subplot(3,nx,3+nx*2)
hold off
for bn = 1:3
    eval(['data = mmburst_' num2str(bn) ';']);
    plot([0 circ_diam2], [data(20) data(1:10)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmburst_all;
plot([0 circ_diam2], [data(20) data(1:10)], 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 50 0 ylim2])
xlabel('Centre diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Centre: total response')

subplot(3,nx,4)
hold off
for bn = 1:3
    eval(['data = mmtrans_' num2str(bn) ';']);
    plot(circ_diam2(10)-mask_diam2, data(11:20), 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmtrans_all;
plot(circ_diam2(10)-mask_diam2, data(11:20), 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 42 0 ylim2])
xlabel('Centre diameter - Centre mask diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Surround: transient response')

subplot(3,nx,4+nx)
hold off
for bn = 1:3
    eval(['data = mmstead_' num2str(bn) ';']);
    plot(circ_diam2(10)-mask_diam2, data(11:20), 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmstead_all;
plot(circ_diam2(10)-mask_diam2, data(11:20), 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 42 0 ylim2])
xlabel('Centre diameter - Centre mask diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Surround: steady state response')

subplot(3,nx,4+nx*2)
hold off
for bn = 1:3
    eval(['data = mmburst_' num2str(bn) ';']);
    plot(circ_diam2(10)-mask_diam2, data(11:20), 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(50,bn*15, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(50,bn*15, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmburst_all;
plot(circ_diam2(10)-mask_diam2, data(11:20), 'ko-', 'markerfacecolor', 'k')
text(50,65, ['All   '; 'spikes'], 'color', 'k')

box off
axis([0 42 0 ylim2])
xlabel('Centre diameter - Centre mask diameter (^o)')
ylabel('Rate (Hz)')
text(5, ylim2, 'Surround: total response')













