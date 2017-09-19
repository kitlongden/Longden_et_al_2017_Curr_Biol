
direction_flag = 1;

h_fig4 = figure(4);
set(h_fig4,'color','w','Position',[0 200 930 750]) % 
nx = 7;
ny = 4;

N_cells = 10;
tf = [0 1 2 4 7 10 13 17 25];
N_tf = length(tf);
load ../'05 Motion in stimulus center and surround'/r_burst_data_VT1.mat r_burst_1 r_burst_2 r_burst_3
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

colormap(macaw_grey)
% colormap(lin_grey_map)
% colormap default

%% Plot colormaps

ClimHi1 = 65;
ClimHi2 = 80;

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
% colorbar
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
% colorbar
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

intf = 5;
data = squeeze(mr_burst1(:,:,intf)+mr_burst2(:,:,intf));
mn_b12 = mean(data);
se_b12 = std(data)/sqrt(N_cells);
intf = 3;
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
axis([-0.1 5.1 0 85])
set(gca,'xtick',tf2,'xticklabel',{'0' '1' '2' '4' '7' '10' '' '17' '25'})
xlabel('Stimulus surround velocity (cycles/s)')
set(gca,'ytick',0:10:90,'yticklabel',{'0' '' '20' '' '40' '' '60' '' '80' ''})
ylabel('Spike rate (Hz)')
set(gca,'tickdir','out')
text(1,75,{'2 cycles/s motion in' 'stimulus center'})
text(-1,90,'F','fontweight','bold','fontsize',14)

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
axis([-0.1 5.1 0 85])
set(gca,'xtick',tf2,'xticklabel',{'0' '1' '2' '4' '7' '10' '' '17' '25'})
xlabel('Stimulus center velocity (cycles/s)')
set(gca,'ytick',0:10:90,'yticklabel',{'0' '' '20' '' '40' '' '60' '' '80' ''})
ylabel('Spike rate (Hz)')
set(gca,'tickdir','out')

legend('3+ spike bursts', ['Single spikes & 2 spike bursts'], 'location', 'southwest')
legend boxoff
text(1,75,{'7 cycles/s motion in' 'stimulus surround'})
text(-1,90,'C','fontweight','bold','fontsize',14)

%% Plot stimulus diagram for centre
%subplot(ny,nx,[1 2 8 9])
% subplot(5,9,[1 2 10 11])
subplot(6,9,[1 2 10 11])
hold off
lgc = 0.75*[1 1 1];

% Plot map
load ../'10 Figure 1'/map_workspace1.mat ix2 x2 iy2 y2 x2 x3 x4 x5 y2 y3 y4 y5
line([ix2 x2]', [iy2 y2]','color',lgc,'linewidth',2.0)
hold on
% plot(ix2,iy2,'o','color',lgc,'markerfacecolor',lgc,'markersize',4);
fill([x2 x3 x4 x5]',[y2 y3 y4 y5]', lgc, 'edgecolor', lgc)
clear ix2 x2 iy2 y2 x2 x3 x4 x5 y2 y3 y4 y5
load ../'10 Figure 1'/map_workspace2.mat x x2 y y2 x3 x4 x5 y3 y4 y5
line([x x2]', [y y2]','color',lgc,'linewidth',2.0)
fill([x3 x4 x5]',[y3 y4 y5]', lgc, 'edgecolor', lgc)
clear x x2 y y2 x3 x4 x5 y3 y4 y5
load ../'10 Figure 1'/map_workspace3.mat x_RF x_RF2 y_RF y_RF2 x_RF3 x_RF4 x_RF5 y_RF3 y_RF4 y_RF5
line([x_RF x_RF2]', [y_RF y_RF2]','color',lgc,'linewidth',2.0)
fill([x_RF3 x_RF4 x_RF5]',[y_RF3 y_RF4 y_RF5]', lgc, 'edgecolor', lgc)
clear x_RF x_RF2 y_RF y_RF2 x_RF3 x_RF4 x_RF5 y_RF3 y_RF4 y_RF5

viscircles([60 -45],24,'color','k','linewidth',1)
hold on
viscircles([60 -45],12,'color','k','linewidth',1)

axis equal
box off
set(gca,'tickdir', 'out')
set(gca,'ytick',-90:15:-15,'yticklabel',{'' '-75^o' '' '-45^o' '' '-15^o' ''})
set(gca,'xtick',0:15:120,'xticklabel',{'0^o' '' '30^o' '' '60^o' '' '90^o' '' '120^o'})
axis([15 105 -90 -15])
% axis equal
xlabel('Azimuth')
ylabel('Elevation')
text(30,5,{'Local preferred direction'},'color',lgc)
text(-8,3,'A','fontweight', 'bold', 'fontsize',14)
text(51,-45,{'Stimulus' 'center' },'fontsize',7)
text(51,-61.5,{'Stimulus' 'surround' },'fontsize',7)

%% Plot direction tuning

if (direction_flag == 1)
    load ../'06 Direction tuning'/r_burst_data.mat r_burst* r_spont*
    N_cells = 6;
    % Mean response spike rates
    mburst_3  = mean(r_burst_3,3);
    mburst_12 = mean(r_burst_1+r_burst_2 ,3);
    % Spontaneous spike rates
    mspont_3  = mean(r_spont_3 ,3);
    mspont_12 = mean(r_spont_1+r_spont_2 ,3);
    % Normalise
    for nc = 1:N_cells
        mburst_3(nc,:)  = mburst_3(nc,:)  - mean(mspont_3(nc,:));
        mburst_12(nc,:) = mburst_12(nc,:) - mean(mspont_12(nc,:));
        mburst_3(nc,:)  = mburst_3(nc,:)/max(mburst_3(nc,:));
        mburst_12(nc,:)  = mburst_12(nc,:)/max(mburst_12(nc,:));
    end
    % Mean and SE
    mmburst_3 = mean(mburst_3);
    semburst_3 = std(mburst_3)/sqrt(N_cells);
    mmburst_12 = mean(mburst_12);
    semburst_12 = std(mburst_12)/sqrt(N_cells);
    
    % Plotting variables
    d_angles = 0:22.5:359;
    N_directions = length(d_angles);
    bucol = zeros(3,3);
    bucol(1,:) = [0 0 1];
    bucol(2,:) = [0.5 0 0.5];
    bucol(3,:) = [1 0 0];
    mks = 2;
    
    subplot(6,9,[28 29 37 38])
    hold off
    plot([d_angles 360], [mmburst_3(1:N_directions) mmburst_3(1)], '^-', 'color', bucol(3,:), 'markerfacecolor', bucol(3,:),'markersize',mks,'linewidth',1.0)
    hold on
    plot([d_angles 360], [mmburst_12(1:N_directions) mmburst_12(1)], 'o-', 'color', bucol(1,:), 'markerfacecolor', bucol(1,:),'markersize',mks,'linewidth',1.0)
     errorbar([d_angles 360], [mmburst_3(1:N_directions) mmburst_3(1)], [semburst_3(1:N_directions) semburst_3(1)]/sqrt(N_cells),'^-', 'color', bucol(3,:), 'markerfacecolor', bucol(3,:),'markersize',mks,'linewidth',0.5)
    errorbar([d_angles 360], [mmburst_12(1:N_directions) mmburst_12(1)], [semburst_12(1:N_directions) semburst_12(1)]/sqrt(N_cells),'o-', 'color', bucol(1,:), 'markerfacecolor', bucol(1,:),'markersize',mks,'linewidth',0.5)
    
    box off
    axis([0 360 -0.3 1])
    set(gca,'tickdir','out')
    set(gca,'xtick', 0:90:360,'xticklabel',{'0^o', '90^o', '180^o', '270^o', '360^o'})
    set(gca,'ytick',-0.3:0.1:1,'yticklabel',{'' '-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6' '' '0.8' '' '1.0'})
    ylabel({'Normalised (Spike rate - ' 'spontaneous spike rate) (Hz)'})
    xlabel('Motion direction')
    
    text(-90,1.1,'D','fontweight','bold','fontsize',14)
    legend(['3+ spike bursts'], ['Single spikes & ' 10 '2 spike bursts'], 'location', 'northwest')
    legend boxoff
    
    clear N_cells
end

