%% Figure 2
%
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

plot_flag = 0;
calc_flag = 0;
all_cells_flag = 0;

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


% MEAN RESPONSES FOR ALL SPATIAL & TEMPORAL FREQUENCIES & EVERY CELL
t = (1:N_samples)/samp_freq;
ndcol = zeros(N_diam,3);
ndcol(:,1) = (1:N_diam)/N_diam;
ndcol(:,3) = (N_diam:-1:1)/N_diam;
lgc = 0.5*[1 1 1];
lgc2 = 0.75*[1 1 1];

bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];

load ../'02 Size tuning'/r_burst_data.mat afbarray_*
h = figure(2);
set(h,'Position', [60, 200, 1200, 400])
set(h,'color',[1 1 1])
nx = 4;
dy = 150;
dyt = num2str(dy/2);
ylim = 1550;

%%  PSTHs
ny = 5;
nx = 4;
for nd = [4 10]
    if nd == 4
        subplot(ny,nx,[10 18])
    else
        subplot(ny,nx,[11 19])
    end
    hold off
    % 1/2 spike bursts
    bn = 3;
    eval(['data3 = squeeze(afbarray_' num2str(bn) '(nd,:));'])
    plot(t, data3, 'color', bucol(bn,:),'linewidth',1)
    hold on
    bn = 2;
    eval(['data2 = squeeze(afbarray_' num2str(bn) '(nd,:));'])
    plot(t, data2, 'color', bucol(bn,:),'linewidth',1)
    % 3+ spike bursts
    bn = 1;
    eval(['data1 = squeeze(afbarray_' num2str(bn) '(nd,:));'])
    plot(t, data1, 'color', bucol(bn,:),'linewidth',1)
    % All spikes
    if (nd == 4)
        legend('3+ spike bursts', '2 spike bursts', 'Single spikes', 'location', 'northeast')
        legend boxoff
    end
    
    axis([0 1.5 0 105])
    box off
    set(gca,'tickdir','out')
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    set(gca,'xtick',[0 0.5 1 1.5])
    set(gca,'ytick',0:20:140)
    
    if (nd == 4)
        text(0.05, 105, ['24^o stimulus diameter' 10 'N=10 cells'])
        text(-0.15, 110, 'B2','fontweight','bold','fontsize',14)
    elseif (nd == 10)
        text(0.05, 105, ['48^o stimulus diameter' 10 'N=10 cells'])
        text(-0.15, 110, 'C2','fontweight','bold','fontsize',14)
    end
end

%% SIZE TUNING - circles
load ../'02 Size tuning'/r_burst_data.mat r_trans* r_stead* r_burst*
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


ylim2 = 110;



% Plot centre errorbar
subplot(1,4,1)
mk = ['o', 's', '^', 'd'];
hold off
data = mburst_all;
plot([0 circ_diam2], [mean(data(:,20)) mean(data(:,1:10))], 'k-', 'marker', mk(4), 'markerfacecolor', 'k','linewidth',1.5)
hold on
for bn = 3:-1:1
    eval(['data = mburst_' num2str(bn) ';']);
    plot([0 circ_diam2], [mean(data(:,20)) mean(data(:,1:10))], '-', 'marker', mk(bn), 'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:),'linewidth',1.5)
    hold on
end
legend('All spikes', '3+ spike bursts', '2 spike bursts', 'Single spikes', 'location', 'southeast')
legend boxoff
for bn = 1:3
    eval(['data = mburst_' num2str(bn) ';']);
    errorbar([0 circ_diam2], [mean(data(:,20)) mean(data(:,1:10))], [std(data(:,20)) std(data(:,1:10))]/sqrt(N_cells),'-', 'marker',mk(bn),'color', bucol(bn,:), 'markerfacecolor', bucol(bn,:),'linewidth',1.0)
    hold on
end
data = mburst_all;
errorbar([0 circ_diam2], [mean(data(:,20)) mean(data(:,1:10))], [std(data(:,20)) std(data(:,1:10))]/sqrt(N_cells),'k-','marker',mk(4), 'markerfacecolor', 'k','linewidth',1.0)

box off
set(gca,'tickdir','out')
% axis([0 50 0 ylim2])
axis([0 50 0 90])
xlabel('Stimulus diameter')
ylabel('Spike rate (Hz)')
set(gca,'xtick',[0:12:48],'xticklabel',{'0^o', '12^o', '24^o', '36^o', '48^o'})
text(-10, 90, 'A','fontweight','bold','fontsize',14)

%% Raster plot

raster_flag = 1;
if (raster_flag == 1)
    % N.b. bindata is size:= 20 x 30000, and there are 30 trials:= trial_X.mat
    datdir{1} = '''../Data/Size tuning/14 Apr 2014-2/bindata/''';
    datdir{2} = '''../Data/Size tuning/17 Mar 2014-2/bindata/''';
    datdir{3} = '''../Data/Size tuning/19 Mar 2014-2/bindata/''';
    datdir{4} = '''../Data/Size tuning/24 Mar 2014-2/bindata/''';
    datdir{5} = '''../Data/Size tuning/24 Mar 2014-3/bindata/''';
    datdir{6} = '''../Data/Size tuning/26 Mar 2014-3/bindata/''';
    datdir{7} = '''../Data/Size tuning/19 Mar 2014-6/bindata/''';
    datdir{8} = '''../Data/Size tuning/21 Mar 2014-2/bindata/''';
    datdir{9} = '''../Data/Size tuning/25 Mar 2014-1/bindata/''';
    datdir{10} = '''../Data/Size tuning/27 Mar 2014-2/bindata/''';
    
    N_cells = 10;
    ny = 5;
    nx = 4;
    subplot(ny,nx,[2 6])
    hold off
    % total_activity = zeros(1,N_cells);
    nstim = 4;
    for nc = 9 %:
        % Load spike times
        id = cell(N_size, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(size, number of samples).
            for ns = 1:N_size
                id{ns,nt} = (find(data(ns,:) > 0))';
                id{ns,nt} = id{ns,nt} / 20; % conversion to ms
            end
            % total_activity(nc) = total_activity(nc) + sum(sum(data));
        end
        % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_size_tuning(id, N_size, N_trials, N_samples);
        
        for nt = 1:N_trials
            data = squeeze(burst_array(nstim,nt,:));
            for nb = 1:3
                if (nb < 3)
                    id =  find(data == nb)/20; % ms
                else
                    id = find(data >= nb)/20;
                end
                if (isempty(id) == 0)
                    Nid = length(id);
                    for nid = 1:Nid
                        line([id(nid) id(nid)], [(nt-1) nt]+((4-nb)-1)*(N_trials+20), 'color', bucol(nb,:))
                        hold on
                    end
                end
            end
        end
    end
    axis([0 1500 -20 150])
    set(gca,'tickdir','out')
    set(gca,'xtick',0:500:1500,'xticklabel',{'0', '0.5', '1.0', '1.5'})
    set(gca,'ytick',[0,30, 50,80, 100,130], 'yticklabel',{'0','30', '0','30', '0','30'})
    xlabel('Time (s)')
    ylabel('Trial #')
    text(30,150,['24^o stimulus diameter' 10 'Cell #9'])
    text(1000,140,'Single spikes', 'color', bucol(1,:))
    text(1000,90,'2 spike bursts', 'color', bucol(2,:))
    text(1000,40,'3+ spike bursts', 'color', bucol(3,:))
    axis off
    text(-150, 150, 'B1','fontweight','bold','fontsize',14)
    
    subplot(ny,nx,[3 7])
    hold off
    % total_activity = zeros(1,N_cells);
    nstim = 10;
    for nc = 9%:N_cells
        % Load spike times
        id = cell(N_size, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(size, number of samples).
            for ns = 1:N_size
                id{ns,nt} = (find(data(ns,:) > 0))';
                id{ns,nt} = id{ns,nt} / 20; % conversion to ms
            end
            % total_activity(nc) = total_activity(nc) + sum(sum(data));
        end
        % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_size_tuning(id, N_size, N_trials, N_samples);
        
        for nt = 1:N_trials
            data = squeeze(burst_array(nstim,nt,:));
            for nb = 1:3
                if (nb < 3)
                    id =  find(data == nb)/20; % ms
                else
                    id = find(data >= nb)/20;
                end
                if (isempty(id) == 0)
                    Nid = length(id);
                    for nid = 1:Nid
                        line([id(nid) id(nid)], [(nt-1) nt]+((4-nb)-1)*(N_trials+20), 'color', bucol(nb,:))
                        hold on
                    end
                end
            end
        end
        axis([0 1500 -10 150])
        set(gca,'tickdir','out')
        set(gca,'xtick',0:500:1500,'xticklabel',{'0', '0.5', '1.0', '1.5'})
        set(gca,'ytick',[0,30, 50,80, 100,130],'yticklabel',{'0','30','0', '30','0','30'})
        xlabel('Time (s)')
        ylabel('Trial #')
        text(30,150,['48^o stimulus diameter' 10 'Cell #9'])
        axis off
        text(-150, 150, 'C1','fontweight','bold','fontsize',14)
    end
end

%% Spatial frequency tuning

N_cells = 10;
N_trials = 25;
N_tf = 8;
N_sw = 6;
N_samples = 20000;

tf = [0.1 0.5 2 4 7 10 17 25];
sw = [5.6160   11.1262   13.8113   16.4361   23.8694   30.5406];
sw2 = round(sw);

load ../'03 Spatial wavelength and temporal frequency tuning'/r_burst_data.mat

mm123 = zeros(N_tf,N_sw);
sem123 = zeros(N_tf,N_sw);
mm1 = zeros(N_tf,N_sw);
sem1 = zeros(N_tf,N_sw);
mm2 = zeros(N_tf,N_sw);
sem2 = zeros(N_tf,N_sw);
mm3 = zeros(N_tf,N_sw);
sem3 = zeros(N_tf,N_sw);
for ntf = 1:N_tf    
    % Calculate spatial frequency tuning
    mrburst_123 = mean(r_burst_1(:,:,ntf,:)+r_burst_2(:,:,ntf,:)+r_burst_3(:,:,ntf,:),4);
    mrburst_1 = mean(r_burst_1(:,:,ntf,:),4);
    mrburst_2 = mean(r_burst_2(:,:,ntf,:),4);
    mrburst_3 = mean(r_burst_3(:,:,ntf,:),4);
    mm123(ntf,:) = nanmean(mrburst_123);
    sem123(ntf,:) = nanstd(mrburst_123)/sqrt(N_cells);
    mm1(ntf,:) = nanmean(mrburst_1);
    sem1(ntf,:) = nanstd(mrburst_1)/sqrt(N_cells);
    mm2(ntf,:) = nanmean(mrburst_2);
    sem2(ntf,:) = nanstd(mrburst_2)/sqrt(N_cells);
    mm3(ntf,:) = nanmean(mrburst_3);
    sem3(ntf,:) = nanstd(mrburst_3)/sqrt(N_cells);
end

subplot(1,4,4)

% Broken down by spike burst number
%
hold off
plot(sw,mm123(5,:), 'd-','color','k','linewidth',1.5,'markerfacecolor','k')
hold on
plot(sw,mm3(5,:), '^-','color',bucol(3,:),'linewidth',1.5,'markerfacecolor',bucol(3,:))
plot(sw,mm2(5,:), 's-','color',bucol(2,:),'linewidth',1.5,'markerfacecolor',bucol(2,:))
plot(sw,mm1(5,:), 'o-','color',bucol(1,:),'linewidth',1.5,'markerfacecolor',bucol(1,:))
errorbar(sw,mm123(5,:),sem123(5,:), 'd-','color','k','linewidth',1.5,'markerfacecolor','k')
errorbar(sw,mm3(5,:),sem3(5,:), '^-','color',bucol(3,:),'linewidth',1.5,'markerfacecolor',bucol(3,:))
errorbar(sw,mm2(5,:),sem2(5,:), 's-','color',bucol(2,:),'linewidth',1.5,'markerfacecolor',bucol(2,:))
errorbar(sw,mm1(5,:),sem1(5,:), 'o-','color',bucol(1,:),'linewidth',1.5,'markerfacecolor',bucol(1,:))

box off
set(gca,'tickdir','out')
set(gca,'xtick',0:5:30,'xticklabel',{'0^o' '5^o' '' '15^o' '' '' '30^o'})
axis([0 32 0 45])
ylabel('Spike rate (Hz)')
xlabel('Spatial frequency')
text(-6,45,'D','fontweight','bold','fontsize',14)

legend('All spikes 7 cycles/s', '3+ spike bursts', '2 spike bursts','single spikes', 'location','northwest')
legend boxoff




