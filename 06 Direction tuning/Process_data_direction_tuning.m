%% TF tuning by bursts

calc_flag = 1;

N_cells = 6;
N_trials = 20;
N_directions = 16;
N_samples = 20000;
N_ms = 2000;
grey_duration = 0.5; % 0.25 before and 0.25 after
stim_duration = 0.5;
samp_freq = 20000;

% The data format is as follows:
% data(directions, samples).

datdir{1} = '''../Data/Direction tuning/15 Jan 2014-2/bindata/''';
datdir{2} = '''../Data/Direction tuning/21 Jan 2014-1/bindata/''';
datdir{3} = '''../Data/Direction tuning/21 Jan 2014-5/bindata/''';
datdir{4} = '''../Data/Direction tuning/21 Jan 2014-6/bindata/''';
datdir{5} = '''../Data/Direction tuning/22 Jan 2014-2/bindata/''';
datdir{6} = '''../Data/Direction tuning/22 Jan 2014-3/bindata/''';

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
    
    r_burst_1 = zeros(N_cells, N_directions, N_trials);
    r_burst_2 = zeros(N_cells, N_directions, N_trials);
    r_burst_3 = zeros(N_cells, N_directions, N_trials);
    
    r_trans_1 = zeros(N_cells, N_directions, N_trials);
    r_trans_2 = zeros(N_cells, N_directions, N_trials);
    r_trans_3 = zeros(N_cells, N_directions, N_trials);
    
    r_stead_1 = zeros(N_cells, N_directions, N_trials);
    r_stead_2 = zeros(N_cells, N_directions, N_trials);
    r_stead_3 = zeros(N_cells, N_directions, N_trials);
    
    r_spont_1 = zeros(N_cells, N_directions, N_trials);
    r_spont_2 = zeros(N_cells, N_directions, N_trials);
    r_spont_3 = zeros(N_cells, N_directions, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    gw = gausswin(200,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_directions,N_samples);
    afbarray_2 = zeros(N_directions,N_samples);
    afbarray_3 = zeros(N_directions,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_directions, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(size, number of samples).
            for nd = 1:N_directions
                id{nd,nt} = (find(data(nd,:) > 0))';
                id{nd,nt} = id{nd,nt} / 20; % conversion to ms
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_direction_tuning(id, N_directions, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for nd = 1:N_directions
            for nt = 1:N_trials
                % Total response
                r_burst_1(nc,nd,nt) = length(find(burst_array(nd,nt,samp_1:samp_2) == 1))/samp_duration;
                r_burst_2(nc,nd,nt) = length(find(burst_array(nd,nt,samp_1:samp_2) == 2))/samp_duration;
                r_burst_3(nc,nd,nt) = length(find(burst_array(nd,nt,samp_1:samp_2) >= 3))/samp_duration;
                % Transient
                r_trans_1(nc,nd,nt) = length(find(burst_array(nd,nt,trans_1:trans_2) == 1))/trans_duration;
                r_trans_2(nc,nd,nt) = length(find(burst_array(nd,nt,trans_1:trans_2) == 2))/trans_duration;
                r_trans_3(nc,nd,nt) = length(find(burst_array(nd,nt,trans_1:trans_2) >= 3))/trans_duration;
                % Steady state
                r_stead_1(nc,nd,nt) = length(find(burst_array(nd,nt,stead_1:stead_2) == 1))/stead_duration;
                r_stead_2(nc,nd,nt) = length(find(burst_array(nd,nt,stead_1:stead_2) == 2))/stead_duration;
                r_stead_3(nc,nd,nt) = length(find(burst_array(nd,nt,stead_1:stead_2) >= 3))/stead_duration;
                % Spontaneous
                r_spont_1(nc,nd,nt) = (length(find(burst_array(nd,nt,spont_1:spont_2) == 1))+length(find(burst_array(nd,nt,spont_3:spont_4) == 1)))/spont_duration;
                r_spont_2(nc,nd,nt) = (length(find(burst_array(nd,nt,spont_1:spont_2) == 2))+length(find(burst_array(nd,nt,spont_3:spont_4) == 2)))/spont_duration;
                r_spont_3(nc,nd,nt) = (length(find(burst_array(nd,nt,spont_1:spont_2) >= 3))+length(find(burst_array(nd,nt,spont_3:spont_4) >= 3)))/spont_duration;
            
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
        for nd = 1:N_directions
            fbarray_1{nc}(nd,:) = filtfilt(gw,1,barray_1{nc}(nd,:));
            fbarray_2{nc}(nd,:) = filtfilt(gw,1,barray_2{nc}(nd,:));
            fbarray_3{nc}(nd,:) = filtfilt(gw,1,barray_3{nc}(nd,:));
            afbarray_1(nd,:) = afbarray_1(nd,:) + fbarray_1{nc}(nd,:);
            afbarray_2(nd,:) = afbarray_2(nd,:) + fbarray_2{nc}(nd,:);
            afbarray_3(nd,:) = afbarray_3(nd,:) + fbarray_3{nc}(nd,:);
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
ndcol = zeros(N_directions,3);
ndcol(:,1) = (1:N_directions)/N_directions;
ndcol(:,3) = (N_directions:-1:1)/N_directions;
lgc = 0.5*[1 1 1];
lgc2 = 0.75*[1 1 1];

load r_burst_data.mat afbarray_*
h = figure(1);
set(h,'Position', [60, 20, 750, 1000])
set(h,'color',[1 1 1])
nx = 2;
dy = 60;
dyt = num2str(dy/2);
ylim = 950;
d_angles = 0:22.5:359;
% PSTHs
for bn = 1:3
    subplot(3,nx,1 + (bn-1)*nx)
    hold off
    patch([0.25 0.3 0.3 0.25], [0 0 ylim ylim], lgc2, 'edgecolor', lgc2)
    hold on
    for nd = 1:N_directions
        eval(['data = squeeze(afbarray_' num2str(bn) '(nd,:));'])
        plot(t, data + (nd-1)*dy, 'color', ndcol(nd,:))
        text(1.01,(nd-1)*dy, [num2str(d_angles(nd)) '^o'], 'color', ndcol(nd,:))
    end
    line([0.25 0.25], [0 ylim], 'color', lgc)
    line([0.75 0.75], [0 ylim], 'color', lgc)
    axis([0 1.0 0 ylim])
    box off
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    set(gca,'xtick',[0 0.25 0.75 1.0])
    set(gca,'ytick',0:dy/2:1450, 'yticklabel',{'0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt, '0', dyt} )
    if (bn < 3)
        text(0.825, ylim, [num2str(bn) ' spb'], 'fontweight', 'bold')
    else
        text(0.825, ylim, [num2str(bn) '+ spb'], 'fontweight', 'bold')
    end
    text(1.0,ylim,['Direction'], 'color', 'r')
    text(0.25,ylim*1.05,'Transient response', 'color', lgc)
    text(0,ylim*1.12,'PSTHs of mean responses', 'color', 'k')
end


% Direction TUNING FOR TRANSIENT AND STEADY STATE RESPONSES
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
ylim2 = 100;

subplot(3,nx,2)
hold off
for bn = 1:3
    eval(['data = mmtrans_' num2str(bn) ';']);
    plot([d_angles 360], [data(1:N_directions) data(1)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(370,bn*10, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(370,bn*10, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmtrans_all;
plot([d_angles 360], [data(1:N_directions) data(1)], 'ko-', 'markerfacecolor', 'k')
text(370,40, 'All spikes', 'color', 'k')

box off
axis([0 360 0 ylim2])
set(gca,'xtick', 0:90:360)
ylabel('Rate (Hz)')
xlabel('Direction (^o)')
text(15, ylim2, 'Transient response')

subplot(3,nx,2+nx)
hold off
for bn = 1:3
    eval(['data = mmstead_' num2str(bn) ';']);
    plot([d_angles 360], [data(1:N_directions) data(1)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(370,bn*10, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(370,bn*10, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmstead_all;
plot([d_angles 360], [data(1:N_directions) data(1)], 'ko-', 'markerfacecolor', 'k')
text(370,40, 'All spikes', 'color', 'k')

box off
axis([0 360 0 ylim2])
set(gca,'xtick', 0:90:360)
ylabel('Rate (Hz)')
xlabel('Direction (^o)')
text(15, ylim2, 'Steady state response')

subplot(3,nx,2+2*nx)
hold off
for bn = 1:3
    eval(['data = mmburst_' num2str(bn) ';']);
    plot([d_angles 360], [data(1:N_directions) data(1)], 'o-', 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:))
    hold on
    if (bn < 3)
        text(370,bn*10, [num2str(bn) ' spb'], 'color', bucol(bn,:))
    else
        text(370,bn*10, [num2str(bn) '+ spb'], 'color', bucol(bn,:))
    end
end
data = mmburst_all;
plot([d_angles 360], [data(1:N_directions) data(1)], 'ko-', 'markerfacecolor', 'k')
text(370,40, 'All spikes', 'color', 'k')

box off
axis([0 360 0 ylim2])
set(gca,'xtick', 0:90:360)
ylabel('Rate (Hz)')
xlabel('Direction (^o)')
text(15, ylim2, 'Total response')

