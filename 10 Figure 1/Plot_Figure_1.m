function Plot_Figure_1()

%% Plotting map - panel A - the optic flow map

% Load data
datdir{1} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-2/''';
datdir{2} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/08 Jan 2014-1/''';
datdir{3} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-1/''';
datdir{4} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/10 Jan 2014-1/''';
datdir{5} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/14 Jan 2014-1/''';
datdir{6} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/21 Jan 2014-3/''';
datdir{7} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/16 Apr 2008/''';
datdir{8} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/29 Apr 2008/''';
datdir{9} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/21 Mar 2007/''';
datdir{10} = '''../Data/Local motion receptive fields/VT1/Wide field maps/25 Feb 2014-3/''';
N_cells = length(datdir);
s(1:N_cells) = 0.6;
u_all = zeros(66, N_cells);
v_all = zeros(66, N_cells);
for nc = 1:N_cells
    % Load data
    eval(['load ' datdir{nc} 'Map_data.mat']);
    u_all(:,nc) = u;
    v_all(:,nc) = v;
end
mu = mean(u_all,2);
mv = mean(v_all,2);

% Interpolate
imu = zeros(11,17);
imv = zeros(11,17);
ix  = zeros(11,17);
iy  = zeros(11,17);
for n = 1:11
    ix(n,:) = -120:15:120;
end
for m = 1:17
    iy(1:10,m) = 75:-15:-60;
    iy(11,m) = -70;
end
% Rows
imu(1,:) = interp1(x(1:7),mu(1:7),-120:15:120);
imv(1,:) = interp1(x(1:7),mv(1:7),-120:15:120);
imu(3,:) = interp1(x(8:16),mu(8:16),-120:15:120);
imv(3,:) = interp1(x(8:16),mv(8:16),-120:15:120);
imu(5,:) = interp1(x(17:33),mu(17:33),-120:15:120);
imv(5,:) = interp1(x(17:33),mv(17:33),-120:15:120);
imu(7,:) = interp1(x(34:50),mu(34:50),-120:15:120);
imv(7,:) = interp1(x(34:50),mv(34:50),-120:15:120);
imu(9,:) = interp1(x(51:59),mu(51:59),-120:15:120);
imv(9,:) = interp1(x(51:59),mv(51:59),-120:15:120);
imu(11,:) = interp1(x(60:66),mu(60:66),-120:15:120);
imv(11,:) = interp1(x(60:66),mv(60:66),-120:15:120);
% Columns
for m = 1:17
    imu(:,m) = interp1(iy(1:2:11,m), imu(1:2:11,m),[75:-15:-60,-70]);
    imv(:,m) = interp1(iy(1:2:11,m), imv(1:2:11,m),[75:-15:-60,-70]);
end
% Reshape
ix2 = reshape(ix',1,11*17)';
iy2 = reshape(iy',1,11*17)';
imu2 = reshape(imu',1,11*17)';
imv2 = reshape(imv',1,11*17)';

% Add RF MAPS
RFdatdir{1} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/15 Jan 2014-1/''';
RFdatdir{2} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-2/''';
RFdatdir{3} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-4/''';
RFdatdir{4} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-7/''';
RFdatdir{5} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-1/''';
RFdatdir{6} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-4/''';
RFdatdir{7} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/11 Feb 2014-3/''';
RFdatdir{8} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/25 Feb 2014-4/''';
RFdatdir{9} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/27 Feb 2014-1/''';
N_RF_cells = length(RFdatdir);

u_all_RF = zeros(11, N_RF_cells);
v_all_RF = zeros(11, N_RF_cells);

for nc = 1:N_RF_cells
    % Load data
    eval(['load ' RFdatdir{nc} 'Map_data.mat u v x y']);
    u_all_RF(:,nc) = u;
    v_all_RF(:,nc) = v;
    x_RF = x;
    y_RF = y;
end

mu_RF = mean(u_all_RF,2);
mv_RF = mean(v_all_RF,2);

% Reinstate x and y for wide field map data
eval(['load ' datdir{1} 'Map_data.mat x y']);

% Normalise
imu3 = imu2/max(sqrt(imu2.^2 + imv2.^2));
imv3 = imv2/max(sqrt(imu2.^2 + imv2.^2));
% Scaling
ss1 = 0.5;
ss2 = 10; % 9
ss3 = 1.5;
% Nonlinear scaling
mag3 = sqrt(imu3.^2 + imv3.^2);
imu4 = (imu3./mag3).*mag3.^(ss1);
imv4 = (imv3./mag3).*mag3.^(ss1);

ss = 0.15;
lgc = 0.5*[1 1 1];

h_fig1 = figure(1);
set(h_fig1,'color','w','Position',[100 200 1000 400]) 

            
subplot(1,2,1)
hold off
x_rf_box = [-20 145 145 -20 -20];
y_rf_box = [-35 -35 -125 -125 -35];
vlgc= 0.9*[1 1 1];

% Drawing arrow shaft
x2 = ix2 + ss2*imu4;
y2 = iy2 + ss2*imv4;
line([ix2 x2]', [iy2 y2]','color',lgc,'linewidth',1.0)
% Arrow head
arrow_ang = atan2d(imv4,imu4);
arrow_len = sqrt((ss2*imu4).^2 + (ss2*imv4).^2);
arrow_len = arrow_len.^0.5;
dx = ss3*arrow_len.*cosd(arrow_ang+90)/3;
dy = ss3*arrow_len.*sind(arrow_ang+90)/3;
du = ss3*arrow_len.*cosd(arrow_ang);
dv = ss3*arrow_len.*sind(arrow_ang);
x3 = x2 + dx;
y3 = y2 + dy;
x4 = x2 + du;
y4 = y2 + dv;
x5 = x2 - dx;
y5 = y2 - dy;
hold on
fill([x2 x3 x4 x5]',[y2 y3 y4 y5]', lgc, 'edgecolor', lgc)
lgc = 0.5*[1 1 1];

save map_workspace1.mat x2 ix2 y2 iy2 x2 x3 x4 x5 y2 y3 y4 y5 

% REPEAT FOR UNINTERPOLATED DATA
%
% Normalise
mu2 = mu./max(sqrt(mu.^2 + mv.^2));
mv2 = mv./max(sqrt(mu.^2 + mv.^2));
% Nonlinear scaling
mag2 = sqrt(mu2.^2 + mv2.^2);
mu3 = (mu2./mag2).*mag2.^(ss1);
mv3 = (mv2./mag2).*mag2.^(ss1);
% DrawIng arrow shaft
x2 = x + ss2*mu3;
y2 = y + ss2*mv3;
line([x x2]', [y y2]','color','k','linewidth',1.5)
% Arrow head
ss3 = 02;
arrow_ang = atan2d(mv,mu);
arrow_len = sqrt((ss2*mu3).^2 + (ss2*mv3).^2);
arrow_len = arrow_len.^0.5;
dx = ss3*arrow_len.*cosd(arrow_ang+90)/3;
dy = ss3*arrow_len.*sind(arrow_ang+90)/3;
du = ss3*arrow_len.*cosd(arrow_ang);
dv = ss3*arrow_len.*sind(arrow_ang);
x3 = x2 + dx;
y3 = y2 + dy;
x4 = x2 + du;
y4 = y2 + dv;
x5 = x2 - dx;
y5 = y2 - dy;
hold on
fill([x3 x4 x5]',[y3 y4 y5]', 'k')

save map_workspace2.mat x x2 y y2 x3 y3 x4 y4 x5 y5 



% Repeat for uninterpolated RF data
%
% Normalise
mu_RF2 = mu_RF./max(sqrt(mu_RF.^2 + mv_RF.^2));
mv_RF2 = mv_RF./max(sqrt(mu_RF.^2 + mv_RF.^2));
% Nonlinear scaling
mag_RF2 = sqrt(mu_RF2.^2 + mv_RF2.^2);
mu_RF3 = (mu_RF2./mag_RF2).*mag_RF2.^(ss1);
mv_RF3 = (mv_RF2./mag_RF2).*mag_RF2.^(ss1);
% DrawIng arrow shaft
x_RF2 = x_RF + ss2*mu_RF3;
y_RF2 = y_RF + ss2*mv_RF3;
line([x_RF x_RF2]', [y_RF y_RF2]','color','k','linewidth',1.5)
% Arrow head
ss3 = 02;
arrow_ang_RF = atan2d(mv_RF,mu_RF);
arrow_len_RF = sqrt((ss2*mu_RF3).^2 + (ss2*mv_RF3).^2);
arrow_len_RF = arrow_len_RF.^0.5;
dx_RF = ss3*arrow_len_RF.*cosd(arrow_ang_RF+90)/3;
dy_RF = ss3*arrow_len_RF.*sind(arrow_ang_RF+90)/3;
du_RF = ss3*arrow_len_RF.*cosd(arrow_ang_RF);
dv_RF = ss3*arrow_len_RF.*sind(arrow_ang_RF);
x_RF3 = x_RF2 + dx_RF;
y_RF3 = y_RF2 + dy_RF;
x_RF4 = x_RF2 + du_RF;
y_RF4 = y_RF2 + dv_RF;
x_RF5 = x_RF2 - dx_RF;
y_RF5 = y_RF2 - dy_RF;
hold on
fill([x_RF3 x_RF4 x_RF5]',[y_RF3 y_RF4 y_RF5]', 'k')

save map_workspace3.mat x_RF x_RF2 y_RF y_RF2 x_RF3 x_RF4 x_RF5 y_RF3 y_RF4 y_RF5

% line([-150 150],[0 0],'color',lgc,'linewidth',1.5,'linestyle','--');
set(gca,'tickdir','out')
set(gca,'xtick',-120:30:120, 'xticklabel', {'-120^o' '-90^o' '-60^o' '-30^o' '0^o' '30^o' '60^o' '90^o' '120^o'});
set(gca,'ytick',[-70 -45 -15 0 15 45 75],'yticklabel',{'-70^o' '-45^o' '-15^o' '0^o' '15^o' '45^o' '75^o'});
xlabel(texlabel('Azimuth'));
ylabel(texlabel('Elevation'));
axis equal
axis([-135 135 -90 85])
box on

text(-165, 98, 'A','fontweight','bold','fontsize',14)
text(-135, 98, 'VT1','fontweight','normal','fontsize',14)

%% Plot rotation and translation fits

% h_fig101 = figure(101);
% set(h_fig101,'color','w','Position',[100 200 1000 400]) 

load ../Data/'Local motion receptive fields'/VT1/ipsilateral_optic_flow_fits_calculated.mat
Ncells = 10;

subplot(6,6,[10 16])
hold off
errorbar([-10:0 1 2],[NaN(1,11) mean(rscore) mean(tscore)],[NaN(1,11) std(rscore) std(tscore)]/sqrt(Ncells),'ko','linewidth',1,'markerfacecolor','k')
hold on
plot([1 2],[rscore tscore],'.-','color',lgc)
set(gca,'xtick',[1 2],'xticklabel',{'Rotation' 'Translation'})
set(gca,'ytick',0:0.2:1)
ylabel(['Best fit to' 10 'optic flow'])
box off
axis([0.5 2.5 0 0.75])
set(gca,'ytick', 0:0.25:0.75)
set(gca,'tickdir', 'out')
text(-0.35, 0.75, 'B','fontweight','bold','fontsize',14)


%% Spike trace

load spike_trace_data.mat
t = (1:20000)/20000;


subplot(6,12,[22 23])
hold off
plot(t, trial_1_tf_1_1_data, 'k')
hold on
% line([0 0.1],[-1.4 -1.4],'color','k','linewidth', 2) 
% line([0 0],[-1.4 -0.4],'color','k','linewidth', 2) 
box off
% axis([0.2 0.45 -1.4 1.2])
axis([0 1 -1.4 1.2])
axis off
text(-0.25, 1.4, 'C','fontweight','bold','fontsize',14)

subplot(6,12,34)
tt = 10000;
hold off
plot(t(tt:tt+1600), trial_1_tf_1_1_data(tt:tt+1600), 'k')
hold on
box off
axis([(tt/20000)-0.02 (tt+1600)/20000 -1.4 1.2])
axis off

subplot(6,12,35)
tt = 16500;
hold off
plot(t(tt:tt+1600), trial_1_tf_1_1_data(tt:tt+1600), 'k')
hold on
box off
axis([(tt/20000)-0.02 (tt+1600)/20000 -1.4 1.2])
axis off

%% Plotting panels D and E 
% N.b. more cells were included in the data anlysis for the paper than are
% included in the data set. They are not included here because they were
% not used for any other results. The impact on Figure 1D and E is very
% small. If you would like the extra data please email.

x = cell(98,1);
bn = zeros(98,100);

Num_cells = [6 10 10 10 11 9 10 6 7 10 9];
Num_data_sets = 11;

cell_count = 1;
for n = 1:Num_data_sets
    eval(['load ISI_data_' num2str(n) '.mat isi_pre burst_num'])
    for nc = 1:Num_cells(n)
        x{cell_count} = log10(isi_pre{nc});
        bn(cell_count,:) = burst_num{nc};
        cell_count = cell_count + 1;
    end
end

bin_max = 3.6;
db = 0.05;
bins = 0:db:bin_max;
bins2 = db/2:db:bin_max-db/2;
N_bins = length(bins) - 1;
isi = zeros(N_bins,98);

for nc = 1:98
    N = length(x{nc});
    for nx = 1:N
        log_isi1 = ceil(x{nc}(nx)/db);
        isi(log_isi1,nc) = isi(log_isi1,nc) + 1;
    end
    isi(:,nc) = isi(:,nc) / sum(isi(:,nc));
    bn(nc,:) = bn(nc,:) / sum(bn(nc,:));
end

subplot(6,6,[22 28])
hold off
plot(bins2,mean(isi,2),'k')
hold on
plot(bins2,mean(isi,2)+std(isi,[],2)/sqrt(98),'color',lgc)
plot(bins2,mean(isi,2)-std(isi,[],2)/sqrt(98),'color', lgc)
line([log10(5) log10(5)],[0 0.07],'color','k','linestyle','--')

box off
set(gca,'tickdir','out')
set(gca,'xtick', log10([2:10, 20:10:100, 200:100:1000]))
set(gca,'xticklabel',{'2','','','','','','','','10','','','','','','','','','100','','','','','','','','','1000'}); 
set(gca,'ytick', 0:0.01:0.07, 'yticklabel', {'0' '' '' '' '' '' '' '0.07'}) 
xlabel('ISI (ms)')
ylabel('Frequency')
text(-1.25, 0.07, 'D','fontweight','bold','fontsize',14)
axis([0 3 0 0.07])
    
subplot(6,6,[23.5 29.5])
hold off
bar(1:10,nanmean(bn(:,1:10)),1,'facecolor','w')
hold on
errorbar(1:10,nanmean(bn(:,1:10)),nanstd(bn(:,1:10))/sqrt(98),'.','color','k','markerfacecolor','k')
axis([0.5 6.5 0 0.55])
set(gca,'ytick', 0:0.1:0.5, 'yticklabel', {'0' '' '' '' '' '0.5'})
box off
set(gca,'tickdir','out')
xlabel('Spikes per burst')
ylabel('Frequency')
text(-2, 0.575, 'E','fontweight','bold','fontsize',14)


