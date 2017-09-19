

h_fig41 = figure(41);
set(h_fig41,'color','w','Position',[0 200 930 750]) 
nx = 7;
ny = 4;

%% Panels A and D

clear variables
x = cell(6,1);
bn = zeros(6,100);

Num_cells = 6;
Num_data_sets = 1;

cell_count = 1;
load data_V1.mat isi_pre burst_num
for nc = 1:Num_cells
    x{cell_count} = log10(isi_pre{nc});
    bn(cell_count,:) = burst_num{nc};
    cell_count = cell_count + 1;
end

bin_max = 3.6;
db = 0.05;
bins = 0:db:bin_max;
bins2 = db/2:db:bin_max-db/2;
N_bins = length(bins) - 1;
isi = zeros(N_bins,6);

for nc = 1:6
    N = length(x{nc});
    for nx = 1:N
        log_isi1 = ceil(x{nc}(nx)/db);
        isi(log_isi1,nc) = isi(log_isi1,nc) + 1;
    end
    isi(:,nc) = isi(:,nc) / sum(isi(:,nc));
    bn(nc,:) = bn(nc,:) / sum(bn(nc,:));
end

subplot(6,4,[1 5])
hold off
plot(bins2,mean(isi,2),'k')
hold on
plot(bins2,mean(isi,2)+std(isi,[],2)/sqrt(6),'r')
plot(bins2,mean(isi,2)-std(isi,[],2)/sqrt(6),'r')
line([log10(5) log10(5)],[0 0.12],'color','k','linestyle','--')
axis([log10(1.35)*0 3 0 0.12])

box off
set(gca,'tickdir','out')
set(gca,'xtick', log10([2:10, 20:10:100, 200:100:1000]))
set(gca,'xticklabel',{'2','3','','5','','','','','10','','','','','','','','','100','','','','','','','','','1000'}); 
xlabel('ISI_n (ms)')
ylabel('Frequency')
text(-0.75,0.13,'A','fontweight','bold','fontsize',14)

subplot(6,4,[13 17])
hold off
bar(1:10,nanmean(bn(:,1:10)),1,'facecolor','w')
hold on
errorbar(1:10,nanmean(bn(:,1:10)),nanstd(bn(:,1:10))/sqrt(108),'o','color','k','markerfacecolor','k')
axis([0.5 6.5 0 1.0])
box off
set(gca,'tickdir','out')
xlabel('Spikes per burst')
ylabel('Frequency')
text(-0.75,1.07,'D','fontweight','bold','fontsize',14)

%% Panels B, C, E, F


N_cells = 6;
tf = [0 1 2 4 7 10 13 17 25];
N_tf = length(tf);
load ../'05 Motion in stimulus center and surround'/r_burst_data_V1.mat r_burst_1 r_burst_2 r_burst_3
% Data format is data(Outer N_tf, Inner N_tf, N_samples)

bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];
lgc = [0.5 0.5 0.5];

%% Calculation joint distributions

mr_burst1  = zeros(N_cells,N_tf,N_tf);
mr_burst2  = zeros(N_cells,N_tf,N_tf);
mr_burst3  = zeros(N_cells,N_tf,N_tf);
mr_burst12 = zeros(N_cells,N_tf,N_tf);
mr_rat     = zeros(N_cells,N_tf,N_tf);

for nc = 1:N_cells
    for ontf = 1:N_tf
        for intf = 1:N_tf
            mr_burst1(nc,ontf,intf) = mean(r_burst_1(nc,ontf,intf,:));
            mr_burst2(nc,ontf,intf) = mean(r_burst_2(nc,ontf,intf,:));
            mr_burst3(nc,ontf,intf) = mean(r_burst_3(nc,ontf,intf,:));
            % Sum of single spikes, 2 spike bursts
            mr_burst12(nc,ontf,intf) = mean(r_burst_1(nc,ontf,intf,:)+r_burst_2(nc,ontf,intf,:));
            % Ratio
            mr_rat(nc,ontf,intf) = mean(r_burst_3(nc,ontf,intf,:)./(r_burst_1(nc,ontf,intf,:)+r_burst_2(nc,ontf,intf,:)));
        end
    end
end

mb1  = squeeze(mean(mr_burst1));
mb2  = squeeze(mean(mr_burst2));
mb3  = squeeze(mean(mr_burst3));
mb12 = squeeze(mean(mr_burst12));
% mr   = squeeze(nanmean(mr_rat));
mr = mb3./mb12;
% tfi = 0:0.5:25;
% tf2i = 0:0.1:5;
tfi = 0:25;
tf2i = 0:0.2:5;
scale_lim = length(tf2i);
tick_fact = max(tf2i);

[x,y] = meshgrid(tf,tf);
[xi,yi] = meshgrid(tfi,tfi);

mb1i  = interp2(x,y,mb1,xi,yi);
mb2i  = interp2(x,y,mb2,xi,yi);
mb3i  = interp2(x,y,mb3,xi,yi);
mb12i = interp2(x,y,mb12,xi,yi);
mri   = interp2(x,y,mr,xi,yi);

tf2 = sqrt(tf);
[x2,y2] = meshgrid(tf2,tf2);
[x2i,y2i] = meshgrid(tf2i,tf2i);

mb1i2  = interp2(x2,y2,mb1,x2i,y2i);
mb2i2  = interp2(x2,y2,mb2,x2i,y2i);
mb3i2  = interp2(x2,y2,mb3,x2i,y2i);
mb12i2 = interp2(x2,y2,mb12,x2i,y2i);
mri2   = interp2(x2,y2,mr,x2i,y2i);

mb1i  = flipud(mb1i);
mb2i  = flipud(mb2i);
mb3i  = flipud(mb3i);
mb12i = flipud(mb12i);
mri   = flipud(mri);

mb1i2  = flipud(mb1i2);
mb2i2  = flipud(mb2i2);
mb3i2  = flipud(mb3i2);
mb12i2 = flipud(mb12i2);
mri2   = flipud(mri2);

%% Color maps

% Linear grey map with saturation at ends
lin_grey_map = zeros(80,3);
id = 0.2 + 0.8*(1./(1 + exp(-0.15*((1:80) - 40))));
lin_grey_map(:,1) = id;
lin_grey_map(:,2) = id;
lin_grey_map(:,3) = id;

% Macaw map: yellow for >0.5, blue for <0.5, grey for 0.5.
macaw_map = zeros(80,3);
macaw_map(:,1) = ([0:79]/79);
macaw_map(:,2) = ([0:79]/79);
macaw_map(:,3) = (79:-1:0)/79;

% Macaw grey map: yellow for >0.6, blue for <0.4, grey for 0.4-0.6
macaw_grey = zeros(80,3);
macaw_grey(:,1) = (0:79)/79;
macaw_grey(:,2) = (0:79)/79;
% macaw_grey(1:31,3) = ((30:-1:0)+20)/50;
% macaw_grey(32:48,3) = (32:48)/80;
% macaw_grey(49:80,3) = (31:-1:0)/52;
macaw_grey(1:35,3) = ((34:-1:0)+26)/60;
macaw_grey(36:45,3) = (35:44)/79;
macaw_grey(46:80,3) = (34:-1:0)/60;

% colormap(macaw_grey)
colormap(lin_grey_map)
% colormap default

%% Plot colormaps

ClimHi1 = 160;
ClimHi2 = 20;

% subplot(ny,nx,[17 18 24 25])
subplot(6,7,[17 18 24 25]+7)
imagesc(mb12i2)
box off

axis([1 scale_lim 1 scale_lim])
set(gca,'tickdir', 'out')
set(gca,'xtick',1+tf2*tick_fact,'xticklabel',{'0' '1' '2' '4' '7' '10' '13' '17' '25'})
yt = 1+tf2*tick_fact;
yt = cumsum([1 fliplr(diff(yt))]);
set(gca,'ytick',yt,'yticklabel',{'25' '17' '13' '10' '7' '4' '2' '1' '0'})
axis equal
xlabel({'Stimulus center velocity (cycles/s)'})
% ylabel({'Stimulus surround velocity (cycles/s)'})
ylabel({'Stimulus                ' 'surround                ' 'velocity                ' '(cycles/s)                '},'rot',0)
set(gca, 'CLim', [0, ClimHi1]);
line([8 8],[0 26],'color','w','linestyle',':','linewidth',2)
line([0 26],[yt(5) yt(5)],'color','w','linestyle','--','linewidth',1.5)
title('Single spikes & 2 spike bursts')
text(-10,-0.7,'E','fontweight','bold','fontsize',14)

% subplot(ny,nx,[3 4 10 11])
subplot(6,7,[3 4 10 11])
imagesc(mb3i2)
box off
axis([1 scale_lim 1 scale_lim])
set(gca,'tickdir', 'out')
set(gca,'xtick',1+tf2*tick_fact,'xticklabel',{'0' '1' '2' '4' '7' '10' '13' '17' '25'})
set(gca,'ytick',yt,'yticklabel',{'25' '17' '13' '10' '7' '4' '2' '1' '0'})
axis equal
xlabel({'Stimulus center velocity (cycles/s)'})
% ylabel({'Stimulus surround velocity (cycles/s)'})
ylabel({'Stimulus                ' 'surround                ' 'velocity                ' '(cycles/s)                '},'rot',0)
set(gca, 'CLim', [0, ClimHi2]);
line([8 8],[0 26],'color','w','linestyle',':','linewidth',2)
line([0 26],[yt(5) yt(5)],'color','w','linestyle','--','linewidth',1.5)
title('3+ Spike bursts')
text(-10,-0.7,'B','fontweight','bold','fontsize',14)

% subplot(ny,nx,[19 26])
subplot(6,7,[19 26]+7)
set(gca, 'CLim', [0, ClimHi1]);
h_cb = colorbar;
ylabel(h_cb,'Spike rate (Hz)')
axis off

% subplot(ny,nx,[5 12])
subplot(6,7,[5 12])
set(gca, 'CLim', [0, ClimHi2]);
h_cb = colorbar;
ylabel(h_cb,'Spike rate (Hz)')
axis off

%% Plot transects

intf = 3;
data = squeeze(mr_burst1(:,:,intf)+mr_burst2(:,:,intf));
mn_b12 = mean(data);
se_b12 = std(data)/sqrt(N_cells);
data = squeeze(mr_burst3(:,:,intf));
mn_b3 = mean(data);
se_b3 = std(data)/sqrt(N_cells);

% subplot(6,nx,[17 18 24 25]+3+nx)
subplot(6,8,[31 32 39 40])
hold off
plot(tf2,mn_b3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0)
hold on
plot(tf2,mn_b12,'o-','color',(bucol(1,:)+bucol(2,:))/2,'markerfacecolor',(bucol(1,:)+bucol(2,:))/2,'linewidth',1.0)
errorbar(tf2,mn_b12,se_b12,'o-','color',(bucol(1,:)+bucol(2,:))/2,'markerfacecolor',(bucol(1,:)+bucol(2,:))/2,'linewidth',1.0)
errorbar(tf2,mn_b3,se_b3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0)

box off
axis([-0.1 5.1 0 125])
set(gca,'xtick',tf2,'xticklabel',{'0' '1' '2' '4' '7' '10' '' '17' '25'})
xlabel('Stimulus surround velocity (cycles/s)')
set(gca,'ytick',0:10:125,'yticklabel',{'0' '' '20' '' '40' '' '60' '' '80' '' '100' '' '120'})
ylabel('Spike rate (Hz)')
set(gca,'tickdir','out')
text(1,75,{'2 cycles/s motion in' 'stimulus center'})
text(-1.25,130,'F','fontweight','bold','fontsize',14)

ontf = 5;
data = squeeze(mr_burst1(:,ontf,:)+mr_burst2(:,ontf,:));
mn_b12 = mean(data);
se_b12 = std(data)/sqrt(N_cells);
data = squeeze(mr_burst3(:,ontf,:));
mn_b3 = mean(data);
se_b3 = std(data)/sqrt(N_cells);

% subplot(ny,nx,[3 4 10 11]+3)
subplot(6,8,[7 8 15 16])
hold off
plot(tf2,mn_b3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0)
hold on
plot(tf2,mn_b12,'o-','color',(bucol(1,:)+bucol(2,:))/2,'markerfacecolor',(bucol(1,:)+bucol(2,:))/2,'linewidth',1.0)
errorbar(tf2,mn_b12,se_b12,'o-','color',(bucol(1,:)+bucol(2,:))/2,'markerfacecolor',(bucol(1,:)+bucol(2,:))/2,'linewidth',1.0)
errorbar(tf2,mn_b3,se_b3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0)

box off
axis([-0.1 5.1 0 125])
set(gca,'xtick',tf2,'xticklabel',{'0' '1' '2' '4' '7' '10' '' '17' '25'})
xlabel('Stimulus center velocity (cycles/s)')
set(gca,'ytick',0:10:125,'yticklabel',{'0' '' '20' '' '40' '' '60' '' '80' '' '100' '' '120'})
ylabel('Spike rate (Hz)')
set(gca,'tickdir','out')

legend('3+ spike bursts', ['Single spikes & 2 spike bursts'], 'location', 'southwest')
legend boxoff
text(1,75,{'7 cycles/s motion in' 'stimulus surround'})
text(-1.25,130,'C','fontweight','bold','fontsize',14)

