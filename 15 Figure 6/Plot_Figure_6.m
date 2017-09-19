function Plot_Figure_6()

h_fig6 = figure(6);
set(h_fig6,'color','w','Position',[100 200 1000 600]) % 800 wide for 17.6 cm (2 columns, J Neuro)

disp('VT3')
datdir{1} = '''../Data/Local motion receptive fields/VT3/Wide field maps/30 Apr 2008-1/''';
datdir{2} = '''../Data/Local motion receptive fields/VT3/Wide field maps/12 Feb 2014-1/''';
datdir{3} = '''../Data/Local motion receptive fields/VT3/Wide field maps/26 Feb 2014-1/''';
datdir{4} = '''../Data/Local motion receptive fields/VT3/Wide field maps/11 Mar 2014-3/''';
datdir{5} = '''../Data/Local motion receptive fields/VT3/Wide field maps/30 Apr 2008-1/''';
datadir = '../''01 Local motion receptive fields''/Data/VT3/Wide'' field maps''/';
fdir{1} = '30'' Apr 2008-1''/';
fdir{2} = '12'' Feb 2014-1''/';
fdir{3} = '26'' Feb 2014-1''/';
fdir{4} = '11'' Mar 2014-3''/';
fdir{5} = '30'' Apr 2008-1''/';
plot_map_it(datdir,2);
clear datdir fdir datadir

disp('VT2')
datdir{1} = '''../Data/Local motion receptive fields/VT2/Wide field maps/04 Jun 2007-1/''';
datdir{2} = '''../Data/Local motion receptive fields/VT2/Wide field maps/11 Dec 2007-1/''';
datdir{3} = '''../Data/Local motion receptive fields/VT2/Wide field maps/26 Feb 2014-3/''';
datadir = '../Data/''Local motion receptive fields''/VT2/Wide'' field maps''/';
fdir{1} = '04'' Jun 2007-1''/';
fdir{2} = '11'' Dec 2007-1''/';
fdir{3} = '26'' Feb 2014-3''/';
plot_map_it(datdir,1);
clear datdir fdir datadir

disp('Hu')
% Ventral visual field
datdir{1} = '''../Data/Local motion receptive fields/Hu/Wide field maps/08 Apr 2008-1/''';
datdir{2} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Feb 2007-2/''';
% Dorsal and ventral
datdir{3} = '''../Data/Local motion receptive fields/Hu/Wide field maps/19 Feb 2007-1/''';
datdir{4} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Nov 2006-1/''';
% Dorsal visual field
datdir{5} = '''../Data/Local motion receptive fields/Hu/Wide field maps/09 Nov 2009-1/''';
datdir{6} = '''../Data/Local motion receptive fields/Hu/Wide field maps/10 Nov 2009-1/''';
datdir{7} = '''../Data/Local motion receptive fields/Hu/Wide field maps/21 Feb 2007-1/''';
datdir{8} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Mar 2007-1/''';
datadir = '../Data/''Local motion receptive fields''/Hu/Wide'' field maps''/';
fdir{1} = '08'' Apr 2008-1''/';
fdir{2} = '28'' feb 2007-2''/';
fdir{3} = '19'' Feb 2007-1''/';
fdir{4} = '28'' Nov 2006-1''/';
fdir{5} = '09'' Nov 2009-1''/';
fdir{6} = '10'' Nov 2009-1''/';
fdir{7} = '21'' Feb 2007-1''/';
fdir{8} = '28'' Mar 2007-1''/';
plot_map_it(datdir,3);
clear datdir fdir datadir

load ../Data/'Local motion receptive fields'/VT3/contralateral_optic_flow_fits_calculated.mat tscore rscore taxes
tax{1} = taxes*pi/180;
tsc{1} = tscore;
rsc{1} = rscore;
load ../Data/'Local motion receptive fields'/VT2/contralateral_optic_flow_fits_calculated.mat tscore rscore taxes
tax{2} = taxes*pi/180;
tsc{2} = tscore;
rsc{2} = rscore;
load ../Data/'Local motion receptive fields'/Hu/contralateral_optic_flow_fits_calculated.mat tscore rscore taxes
tax{3} = taxes*pi/180;
tsc{3} = tscore;
rsc{3} = rscore;
load VT1_ipsilateral_optic_flow_fits_calculated.mat tscore rscore taxes
tax{4} = taxes*pi/180;
tsc{4} = tscore;
rsc{4} = rscore;
lgc = 0.5*[1 1 1];

% Plot fits to rotation and translation optic flow
% subplot(6,2,8)
panid = [22 28 34]+6;
yval = 0.8;
xval = 0.6;
for n = 1:3
    subplot(7,6,panid(n))
    hold off
    plot([1 2],[rsc{n} tsc{n}],'.-','color',lgc)
    hold on
    errorbar([-10:0 1 2],[NaN(1,11) mean(rsc{n}) mean(tsc{n})],[NaN(1,11) std(rsc{n}) std(tsc{n})]/sqrt(length(rsc{n})),'ok-','markerfacecolor','k','linewidth',1)
    box off
    set(gca,'tickdir', 'out')
    axis([0.5 2.25 0 1])
    set(gca,'xtick',[1 2])
    set(gca,'ytick',0:0.25:1.0,'yticklabel',{'0', '', '0.5', '', '1.0'})
    if (n == 1)
        text(xval,yval,'VT2','fontsize',12)
        % set(gca,'xticklabel',{'' ''})
        set(gca,'xticklabel',{'Rot' 'Trans'})
        text(-0.15,1.35,'D','fontweight','bold','fontsize',14)
    elseif (n == 2)
        ylabel(['Best fit to optic flow'])
        text(xval,yval,'VT3','fontsize',12)
        % set(gca,'xticklabel',{'' ''})
        set(gca,'xticklabel',{'Rot' 'Trans'})
    elseif (n == 3)
        text(xval,yval,'Hu','fontsize',12)
        set(gca,'xticklabel',{'Rot' 'Trans'})
    end
end

% Plot axes of fits
sym = ['d' '^' 's' 'o'];
pcol = ['k' 'b' 'g' 'r'];
mks1 = 4;
mks2 = 8;
% subplot(8,4,[27 31])
subplot(2,3,6)
hold off
polarplot(mean(tax{4}(:,1)),mean(cos(tax{4}(:,2))),'marker',sym(4),'color',pcol(4),'markerfacecolor',pcol(4),'linestyle','none','markersize',mks2);
hold on
polarplot(mean(tax{1}(:,1)),mean(cos(tax{1}(:,2))),'marker',sym(1),'color',pcol(1),'markerfacecolor',pcol(1),'linestyle','none','markersize',mks2);
polarplot(mean(tax{2}(:,1)),mean(cos(tax{2}(:,2))),'marker',sym(2),'color',pcol(2),'markerfacecolor',pcol(2),'linestyle','none','markersize',mks2);
polarplot(mean(tax{3}(:,1)),mean(cos(tax{3}(:,2))),'marker',sym(3),'color',pcol(3),'markerfacecolor',pcol(3),'linestyle','none','markersize',mks2);
polarplot(tax{4}(:,1),cos(tax{4}(:,2)),'marker',sym(4),'color',pcol(4),'markerfacecolor','w','linestyle','none','markersize',mks1);
polarplot(tax{1}(:,1),cos(tax{1}(:,2)),'marker',sym(1),'color',pcol(1),'markerfacecolor','w','linestyle','none','markersize',mks1);
polarplot(tax{2}(:,1),cos(tax{2}(:,2)),'marker',sym(2),'color',pcol(2),'markerfacecolor','w','linestyle','none','markersize',mks1);
polarplot(tax{3}(:,1),cos(tax{3}(:,2)),'marker',sym(3),'color',pcol(3),'markerfacecolor','w','linestyle','none','markersize',mks1);
set(gca,'ThetaZeroLocation', 'top')
set(gca,'ThetaDir','clockwise')
set(gca,'ThetaTick',0:30:360,'ThetaTickLabels',{'0^o' '' '60^o' '' '120^o' '' '' '' '-120^o' '' '-60^o' ''})
set(gca,'RLim',[0 1])
set(gca,'RTick',cosd([90 75 60 45 30 15 0]),'RTicklabel',{'' '75^o' '60^o' '' '30^o' '' '0^o'})
set(gca,'Rcolor',0.6*[1 1 1])
text(0.26,1.0,'Elevation','color',0.6*[1 1 1])
legend('VT1','VT2','VT3','Hu', 'location', 'south')
legend boxoff
text(-0.38,1.2,'Azimuth')
polarplot(mean(tax{4}(:,1))*ones(1,2),[0.5 1.5],'color',pcol(4),'linestyle','-','linewidth',1.5);
polarplot(mean(tax{3}(:,1))*ones(1,2),[0.5 1.5],'color',pcol(3),'linestyle','-','linewidth',1.5);
text(-0.8,1.63,'E','fontweight', 'bold','fontsize',14)
% axis equal


end

function plot_map_it(datdir,cell_no)

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

% Remove measured values from interpolated matrices
for n = 1:length(x)
    id = find(ix2 == x(n) & iy2 == y(n));
    imu2(id) = NaN;
    imv2(id) = NaN;
end

% Reinstate x and y for wide field map data
eval(['load ' datdir{1} 'Map_data.mat x y']);

% Normalise
imu3 = imu2/max(sqrt(imu2.^2 + imv2.^2));
imv3 = imv2/max(sqrt(imu2.^2 + imv2.^2));
% Scaling
ss1 = 0.5; % 0.5
ss2 = 8; % 10
ss3 = 1.5;
% Nonlinear scaling
mag3 = sqrt(imu3.^2 + imv3.^2);
imu4 = (imu3./mag3).*mag3.^(ss1);
imv4 = (imv3./mag3).*mag3.^(ss1);
% imu4 = imu3;
% imv4 = imv3;

ss = 0.15;
lgc = 0.5*[1 1 1];

subplot(2,2,cell_no)
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

% line([-150 150],[0 0],'color',lgc,'linewidth',1.5,'linestyle','--');
set(gca,'tickdir','out')
set(gca,'xtick',-120:30:120, 'xticklabel', {'-120^o' '-90^o' '-60^o' '-30^o' '0^o' '30^o' '60^o' '90^o' '120^o'});
set(gca,'ytick',[-70 -45 -15 0 15 45 75],'yticklabel',{'-70^o' '-45^o' '-15^o' '0^o' '15^o' '45^o' '75^o'});
xlabel(texlabel('Azimuth'));
ylabel(texlabel('Elevation'));
axis equal
axis([-135 135 -90 85])
box on

yval = 98;
if (cell_no == 1)
    text(-165, yval, 'A','fontweight','bold','fontsize',14)
    text(-135, yval, 'VT2','fontweight','normal','fontsize',14)
elseif (cell_no == 2)
    text(-165, yval, 'B','fontweight','bold','fontsize',14)
    text(-135, yval, 'VT3','fontweight','normal','fontsize',14)
elseif (cell_no == 3)
    text(-165, yval, 'C','fontweight','bold','fontsize',14)
    text(-135, yval, 'Hu','fontweight','normal','fontsize',14)
end

end

