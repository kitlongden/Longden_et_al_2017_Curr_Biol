
% Figure 3

N_cells = 10;
N_tf = 9;
dt = 3000;

bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];
lgc = [0.5 0.5 0.5];

h_fig4 = figure(4);
set(h_fig4,'color','w','Position',[100 200 1000 600]) 

% Plot stimulus diagram for STA
subplot(16,10,[14 24])
hold off
viscircles([60 -45],24,'color','k','linewidth',0.5)
hold on

box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-90:15:0,'yticklabel',{'' '-75^o' '' '-45^o' '' '-15^o' ''})
set(gca,'xtick',0:15:120,'xticklabel',{'0^o' '' '30^o' '' '60^o' '' '90^o' '' '120^o'})
% axis([0 120 -90 0])
axis([30 90 -75 -15])
axis equal
xlabel('Azimuth')
ylabel('Elevation')
text(30,5,{'White noise velocity' 'along preferred' 'orientation'})
text(-20,15,'A1','fontweight', 'bold', 'fontsize',14)

% Plot stimulus diagram for centre + surround
subplot(16,10,[54 64])
hold off
viscircles([60 -45],24,'color','k','linewidth',0.5)
hold on

box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-90:15:0,'yticklabel',{'' '-75^o' '' '-45^o' '' '-15^o' ''})
set(gca,'xtick',0:15:120,'xticklabel',{'0^o' '' '30^o' '' '60^o' '' '90^o' '' '120^o'})
% axis([0 120 -90 0])
axis([30 90 -75 -15])
axis equal
xlabel('Azimuth')
ylabel('Elevation')
text(30,-5,{'Motion in preferred', 'direction'})
text(-20,10,'B1','fontweight', 'bold', 'fontsize',14)

% Plot stimulus diagram for centre
subplot(16,10,[94 104])
hold off
viscircles([60 -45],24,'color','k','linewidth',0.5)
hold on
viscircles([60 -45],12,'color','k','linewidth',0.5)

box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-90:15:0,'yticklabel',{'' '-75^o' '' '-45^o' '' '-15^o' ''})
set(gca,'xtick',0:15:120,'xticklabel',{'0^o' '' '30^o' '' '60^o' '' '90^o' '' '120^o'})
% axis([0 120 -90 0])
axis([30 90 -75 -15])
axis equal
xlabel('Azimuth')
ylabel('Elevation')
text(30,-5,{'Motion in center'})
text(-20,10,'C1','fontweight', 'bold', 'fontsize',14)

% Plot stimulus diagram for surround
subplot(16,10,[134 144])
hold off
viscircles([60 -45],24,'color','k','linewidth',0.5)
hold on
viscircles([60 -45],12,'color','k','linewidth',0.5)
box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-90:15:0,'yticklabel',{'' '-75^o' '' '-45^o' '' '-15^o' ''})
set(gca,'xtick',0:15:120,'xticklabel',{'0^o' '' '30^o' '' '60^o' '' '90^o' '' '120^o'})
% axis([0 120 -90 0])
axis([30 90 -75 -15])
axis equal
xlabel('Azimuth')
ylabel('Elevation')
text(30,-5,{'Motion in surround'})
text(-20,10,'D1','fontweight', 'bold', 'fontsize',14)

%% Plot STA stimuli and responses
load datan10.mat datan10 frameshift
t1 = 29200*20; % samples
t2 = 30000*20;
t = (1:600000)/20;
tt1 = 29200/5;  % frame displays
tt2 = 30000/5;
tt = (1:6000)*5;
ffact = 25/16.6;
frameshift = frameshift * ffact;
spike1t = 29.452;
spike2t = 29.828;
dtstim = 0.15;
st1 = spike1t*20000-80;
st2 = spike1t*20000+80;
st3 = spike2t*20000-80;
st4 = spike2t*20000+180;

% Spike traces
subplot(8,5,5)
hold off
plot(t(t1:t2)/1000,datan10(t1:t2),'color',lgc)
hold on
plot(t(st1:st2)/1000,datan10(st1:st2),'color',bucol(1,:))
plot(t(st3:st4)/1000,datan10(st3:st4),'color',bucol(2,:))
box off
axis([29.2 30 -2.5 2.5])
set(gca,'tickdir', 'out')
axis off
text(28.89, 2.5, 'A3', 'fontweight', 'bold', 'fontsize', 14);


% Stimuli
subplot(8,5,10)
hold off
patch([spike1t spike1t spike1t-dtstim spike1t-dtstim],[-25 25 25 -25],bucol(1,:),'facealpha',0.25,'edgecolor','none')
hold on
patch([spike2t spike2t spike2t-dtstim spike2t-dtstim],[-25 25 25 -25],bucol(2,:),'facealpha',0.25,'edgecolor','none')
plot(tt(tt1:tt2)/1000,frameshift(10,tt1:tt2),'k')
box off
axis([29.2 30 -25 25])
set(gca,'tickdir', 'out')
set(gca,'xtick',29.2:0.1:30,'xticklabel',{'' '' '' '29.5' '' '' '' '' '30'})
set(gca,'ytick',-25:25:25)
xlabel('Time (s)')
% ylabel(['Stimulus' 10 'velocity' 10 '(cycles/s)'])
text(28.8,0,['Stimulus' 10 'velocity' 10 '(cycles/s)'])

%% Plot STAs
load ../'04 Random motion stimulus'/sta_data mn_mstim_all se_mstim_all


subplot(4,8,[5 6])


hold off
plot((-dt+1:0)/20,-mn_mstim_all(4,:),'color',bucol(3,:));
hold on
plot((-dt+1:0)/20,-mn_mstim_all(2,:),'color',bucol(2,:));
plot((-dt+1:0)/20,-mn_mstim_all(1,:),'color',bucol(1,:));
% I have edited shadedErrorBar to never do line edges.
% Last input (...'1') is for transparency and using OpenGL as the renderer.
% shadedErrorBar((-dt+1:0)/20,-mn_mstim_all(4,:),se_mstim_all(4,:),{'color',bucol(3,:)},1);
% shadedErrorBar((-dt+1:0)/20,-mn_mstim_all(2,:),se_mstim_all(2,:),{'color',bucol(2,:)},1)
% shadedErrorBar((-dt+1:0)/20,-mn_mstim_all(1,:),se_mstim_all(1,:),{'color',bucol(1,:)},1)

axis([-150 0 -0.35 4.2])
box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-0.5:.5:4, 'yticklabel', {'' '0' '' '1.0' '' '2.0' '' '3.0' '' '4.0'})
set(gca,'xtick',-150:25:0, 'xticklabel', {'-150' '' '-100' '' '-50' '' '0'})
xlabel('Time (ms)')
ylabel({'Stimulus' 'velocity' '(cycles/s)'})
text(-140,4.2,['Spike burst-triggered averages'])
line([-140 -130],[3.5 3.5],'color',bucol(3,:),'linewidth',1.0)
line([-140 -130],[3.0 3.0],'color',bucol(2,:),'linewidth',1.0)
line([-140 -130],[2.5 2.5],'color',bucol(1,:),'linewidth',1.0)
text(-127, 3.55, '3+ spike bursts')
text(-127, 3.05, '2 spike bursts')
text(-127, 2.55, 'Single spikes')
% legend boxoff
text(-185, 4.2, 'A2', 'fontweight', 'bold', 'fontsize', 14);

%% Plot surround TF tuning
% Temporal frequencies are:
tf = [0 1 2 4 7 10 13 17 25];
% Data format is data(Outer N_tf, Inner N_tf, N_samples)
% r_burst_N values are in Hz
% The r_burst_3 is for 3+ spike bursts
% The spike sorting has been updated to catch the ISI<1ms bugs.

load ../'05 Motion in stimulus center and surround'/r_burst_data_VT1.mat r_burst_1 r_burst_2 r_burst_3

ny = 4;
nx = 8;
ni3 = [5 6]+8;
ni2 = [13 14]+8;
ni1 = [21 22]+8;

mr_burst1 = zeros(N_cells,N_tf);
mr_burst2 = zeros(N_cells,N_tf);
mr_burst3 = zeros(N_cells,N_tf);
for nc = 1:N_cells
    for ntf = 1:N_tf
            mr_burst1(nc,ntf) = mean(r_burst_1(nc,ntf,1,:));
            mr_burst2(nc,ntf) = mean(r_burst_2(nc,ntf,1,:));
            mr_burst3(nc,ntf) = mean(r_burst_3(nc,ntf,1,:));
    end
end
mnb1 = mean(mr_burst1);
mnb2 = mean(mr_burst2);
mnb3 = mean(mr_burst3);
seb1 = std(mr_burst1)/sqrt(N_cells);
seb2 = std(mr_burst2)/sqrt(N_cells);
seb3 = std(mr_burst3)/sqrt(N_cells);

subplot(ny,nx,ni1)
hold off
plot(tf,mnb3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:))
hold on
plot(tf,mnb2,'s-','color',bucol(2,:),'markerfacecolor',bucol(2,:))
plot(tf,mnb1,'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(tf,mnb1,seb1,'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(tf,mnb2,seb2,'s-','color',bucol(2,:),'markerfacecolor',bucol(2,:))
errorbar(tf,mnb3,seb3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:))

axis([0 25 0 90])
box off
set(gca,'tickdir','out')
xlabel('Stimulus velocity (cycles/s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',tf)
%set(gca,'xticklabel', {'' '' '' '' '' '' '' '' ''})
set(gca,'ytick',0:10:90)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' ''})
text(-5, 90, 'D2','fontweight', 'bold', 'fontsize',14)



%% Plot centre TF tuning
mr_burst1 = zeros(N_cells,N_tf);
mr_burst2 = zeros(N_cells,N_tf);
mr_burst3 = zeros(N_cells,N_tf);
for nc = 1:N_cells
    for ntf = 1:N_tf
            mr_burst1(nc,ntf) = mean(r_burst_1(nc,1,ntf,:));
            mr_burst2(nc,ntf) = mean(r_burst_2(nc,1,ntf,:));
            mr_burst3(nc,ntf) = mean(r_burst_3(nc,1,ntf,:));
    end
end
mnb1 = mean(mr_burst1);
mnb2 = mean(mr_burst2);
mnb3 = mean(mr_burst3);
seb1 = std(mr_burst1)/sqrt(N_cells);
seb2 = std(mr_burst2)/sqrt(N_cells);
seb3 = std(mr_burst3)/sqrt(N_cells);

subplot(ny,nx,ni2)
hold off
plot(tf,mnb3,'color',bucol(3,:))
hold on
plot(tf,mnb2,'color',bucol(2,:))
plot(tf,mnb1,'color',bucol(1,:))
errorbar(tf,mnb1,seb1,'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(tf,mnb2,seb2,'s-','color',bucol(2,:),'markerfacecolor',bucol(2,:))
errorbar(tf,mnb3,seb3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:))

axis([0 25 0 90])
box off
set(gca,'tickdir','out')
xlabel('Stimulus velocity (cycles/s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',tf)
%set(gca,'xticklabel', {'' '' '' '' '' '' '' '' ''})
set(gca,'ytick',0:10:90)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' ''})
text(-5, 90, 'C2','fontweight', 'bold', 'fontsize',14)

%% Plot centre+surround TF tuning
mr_burst1 = zeros(N_cells,N_tf);
mr_burst2 = zeros(N_cells,N_tf);
mr_burst3 = zeros(N_cells,N_tf);
for nc = 1:N_cells
    for ntf = 1:N_tf
            mr_burst1(nc,ntf) = mean(r_burst_1(nc,ntf,ntf,:));
            mr_burst2(nc,ntf) = mean(r_burst_2(nc,ntf,ntf,:));
            mr_burst3(nc,ntf) = mean(r_burst_3(nc,ntf,ntf,:));
    end
end
mnb1 = mean(mr_burst1);
mnb2 = mean(mr_burst2);
mnb3 = mean(mr_burst3);
seb1 = std(mr_burst1)/sqrt(N_cells);
seb2 = std(mr_burst2)/sqrt(N_cells);
seb3 = std(mr_burst3)/sqrt(N_cells);

subplot(ny,nx,ni3)
hold off
plot(tf,mnb3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1)
hold on
plot(tf,mnb2,'s-','color',bucol(2,:),'markerfacecolor',bucol(2,:),'linewidth',1)
plot(tf,mnb1,'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1)
errorbar(tf,mnb1,seb1,'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(tf,mnb2,seb2,'s-','color',bucol(2,:),'markerfacecolor',bucol(2,:))
errorbar(tf,mnb3,seb3,'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:))

axis([0 25 0 90])
box off
set(gca,'tickdir','out')
xlabel('Stimulus velocity (cycles/s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',tf)
set(gca,'ytick',0:10:90)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' ''})
text(-5, 90, 'B2','fontweight', 'bold', 'fontsize',14)

legend('3+ spike bursts', '2 spike bursts', 'Single spikes', 'location', 'northeast')
legend boxoff

%% Plot PSTHs
% load r_burst_data.mat afbarray_1 afbarray_2 afbarray_3
load ../'05 Motion in stimulus center and surround'/r_burst_data_VT1.mat afbarray_1 afbarray_2 afbarray_3
t = (1:20000)/20;
ntf = 5;

ny = 4;
nx = 5;
ni = [15 10 5]+5;

subplot(ny,nx,ni(1))
hold off
data1 = squeeze(afbarray_1(ntf,1,:));
data2 = squeeze(afbarray_2(ntf,1,:));
data3 = squeeze(afbarray_3(ntf,1,:));
plot(t,data1,'color',bucol(1,:))
hold on
plot(t,data2,'color',bucol(2,:))
plot(t,data3,'color',bucol(3,:))

box off
set(gca,'tickdir', 'out')
axis([0 1000 0 120])
xlabel('Time (s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',[0 250 750 1000])
% set(gca,'xticklabel',{'' '' '' ''})
set(gca,'xticklabel',{'0' '0.25' '0.75' '1'})
set(gca,'ytick',0:10:120)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' '' '' '' '120'})
text(50,120,'7 cycles/s stimulus')
text(-400, 120, 'A3','fontweight', 'bold', 'fontsize',14)

subplot(ny,nx,ni(2))
hold off
data1 = squeeze(afbarray_1(1,ntf,:));
data2 = squeeze(afbarray_2(1,ntf,:));
data3 = squeeze(afbarray_3(1,ntf,:));
plot(t,data1,'color',bucol(1,:))
hold on
plot(t,data2,'color',bucol(2,:))
plot(t,data3,'color',bucol(3,:))

box off
set(gca,'tickdir', 'out')
axis([0 1000 0 120])
xlabel('Time (s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',[0 250 750 1000])
% set(gca,'xticklabel',{'' '' '' ''})
set(gca,'xticklabel',{'0' '0.25' '0.75' '1'})
set(gca,'ytick',0:10:120)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' '' '' '' '120'})
text(-400, 120, 'C3','fontweight', 'bold', 'fontsize',14)

subplot(ny,nx,ni(3))
hold off
data1 = squeeze(afbarray_1(ntf,ntf,:));
data2 = squeeze(afbarray_2(ntf,ntf,:));
data3 = squeeze(afbarray_3(ntf,ntf,:));
plot(t,data1,'color',bucol(1,:))
hold on
plot(t,data2,'color',bucol(2,:))
plot(t,data3,'color',bucol(3,:))

box off
set(gca,'tickdir', 'out')
axis([0 1000 0 120])
xlabel('Time (s)')
ylabel('Spike rate (Hz)')
set(gca,'xtick',[0 250 750 1000])
set(gca,'xticklabel',{'0' '0.25' '0.75' '1'})
set(gca,'ytick',0:10:120)
set(gca,'yticklabel',{'0' '' '20' '' '40' '' '' '' '80' '' '' '' '120'})
text(-400, 120, 'B3','fontweight', 'bold', 'fontsize',14)


