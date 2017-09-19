
%% 
addpath mapfit_fns/
h_fig22 = figure(22);
set(h_fig22,'color','w','Position',[50 0 1200 2000]) % 

%% Plot all maps

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
    eval(['load ' datdir{nc} 'Map_data.mat x y u v']);
    u_all(:,nc) = u;
    v_all(:,nc) = v;
end
mu = mean(u_all,2);
mv = mean(v_all,2);

subplot(2,2,1)
s = 7.5;
s2 = s*1.2;
lw = 1.2;
lgc = [0.5 0.5 0.5];
hold off
for nc = 1:10
    u = u_all(:,nc);
    v = v_all(:,nc);
    for n = 1:66
        mag = sqrt(u(n)^2 + v(n)^2);
        line([x(n) x(n)+s*u(n)/mag],[y(n) y(n)+s*v(n)/mag],'color',lgc)
        hold on
    end
end
for n = 1:66
    u = mu;
    v = mv;
    mag = sqrt(u(n)^2 + v(n)^2);
    plot([x(n) x(n)+s2*u(n)/mag],[y(n) y(n)+s2*v(n)/mag],'color','k','linewidth',lw)
    hold on
end

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

subplot(2,2,1)
for nc = 1:9
    u = u_all_RF(:,nc);
    v = v_all_RF(:,nc);
    for n = 1:11
        mag = sqrt(u(n)^2 + v(n)^2);
        line([x(n) x(n)+s*u(n)/mag],[y(n) y(n)+s*v(n)/mag],'color',lgc)
        hold on
    end
end
for n = 1:11
    u = mu_RF;
    v = mv_RF;
    mag = sqrt(u(n)^2 + v(n)^2);
    plot([x(n) x(n)+s2*u(n)/mag],[y(n) y(n)+s2*v(n)/mag],'color','k','linewidth',lw)
    hold on
end

box on
axis([-135 135 -90 90])
set(gca,'tickdir','out')
set(gca,'xtick',-135:15:135)
set(gca,'ytick',-90:15:90)
set(gca,'xticklabel',{'-135','','','-90','','','-45','','','0','','','45','','','90','','','135'})
set(gca,'yticklabel',{'-90','','','-45','','','0','','','45','','','90'})
xlabel('Azimuth')
ylabel('Elevation')
text(-165,100,'A','fontweight','bold','fontsize',14)



%% Plot fits for cells

ncol = 4;
nrow = 4;
lgc = [0.5 0.5 0.5];
ms1 = 4;
ms2 = 8;

% Ipsilateral fits
load ../Data/'Local motion receptive fields'/VT1/ipsilateral_optic_flow_fits_calculated.mat
subplot(nrow,ncol,9)
hold off
for nc = 1:N_cells
    plot([1, 2],[rscore, tscore],'.-','color',lgc)
    hold on
end
plot([1, 2],[mean(rscore) mean(tscore)],'-','color','k')
errorbar(1,mean(rscore),std(rscore)/sqrt(N_cells),'o-','color','k','markerfacecolor','k')
errorbar(2,mean(tscore),std(tscore)/sqrt(N_cells),'^-','color','k','markerfacecolor','k')
[h,p] = ttest(rscore,tscore);
text(1,0.85,['P = ' num2str(p)])

axis([0.5 2.5 0 1])
box off
ylabel('Fit')
set(gca,'xtick',1:2,'xticklabel',{'Rot', 'Trans'})
set(gca,'Tickdir','out')
set(gca,'ytick',0:0.1:1)
set(gca,'yticklabel',{'0','','','','','0.5','','','','','1'})
text(0.1,1.2,'B1','fontweight','bold','fontsize',14)

subplot(nrow,ncol,10)
hold off
plot(taxes(:,1),taxes(:,2),'k^','markerfacecolor','none','markersize',ms1);
hold on
plot(raxes(:,1),raxes(:,2),'ko','markerfacecolor','none','markersize',ms1);
plot(mean(taxes(:,1)),mean(taxes(:,2)),'k^','markerfacecolor','k','markersize',ms2);
plot(mean(raxes(:,1)),mean(raxes(:,2)),'ko','markerfacecolor','k','markersize',ms2);

axis([-180 180 -90 90])
set(gca,'Tickdir','out')
set(gca,'xtick',-180:30:180)
set(gca,'xticklabel',{'-180','','','90','','','0','','','90','','','180'})
set(gca,'ytick',-90:30:90)
set(gca,'yticklabel',{'-90','','','0','','','90'})
xlabel('Azimuth (^o)')
ylabel('Elevation (^o)')
text(-260,120,'B2','fontweight','bold','fontsize',14)

% Whole visual field fits
load whole_field_optic_flow_fits_calculated.mat
subplot(nrow,ncol,13)
hold off
for nc = 1:N_cells
    plot([1, 2],[rscore, tscore],'.-','color',lgc)
    hold on
end
plot([1, 2],[mean(rscore) mean(tscore)],'-','color','k')
errorbar(1,mean(rscore),std(rscore)/sqrt(N_cells),'o-','color','k','markerfacecolor','k')
errorbar(2,mean(tscore),std(tscore)/sqrt(N_cells),'^-','color','k','markerfacecolor','k')
[h,p] = ttest(rscore,tscore);
text(1,0.85,['P = ' num2str(p)])
 
axis([0.5 2.5 0 1])
box off
ylabel('Fit')
set(gca,'xtick',1:2,'xticklabel',{'Rot', 'Trans'})
set(gca,'Tickdir','out')
set(gca,'ytick',0:0.1:1)
set(gca,'yticklabel',{'0','','','','','0.5','','','','','1'})
text(0.1,1.2,'C1','fontweight','bold','fontsize',14)

subplot(nrow,ncol,14)
hold off
plot(taxes(:,1),taxes(:,2),'k^','markerfacecolor','none','markersize',ms1);
hold on
plot(raxes(:,1),raxes(:,2),'ko','markerfacecolor','none','markersize',ms1);
plot(mean(taxes(:,1)),mean(taxes(:,2)),'k^','markerfacecolor','k','markersize',ms2);
plot(mean(raxes(:,1)),mean(raxes(:,2)),'ko','markerfacecolor','k','markersize',ms2);

axis([-180 180 -90 90])
set(gca,'Tickdir','out')
set(gca,'xtick',-180:30:180)
set(gca,'xticklabel',{'-180','','','90','','','0','','','90','','','180'})
set(gca,'ytick',-90:30:90)
set(gca,'yticklabel',{'-90','','','0','','','90'})
xlabel('Azimuth (^o)')
ylabel('Elevation (^o)')
text(-260,120,'C2','fontweight','bold','fontsize',14)

%% Fits of mean responses of all cells

load Figure_1_mean_data.mat cmap cmap2 camp sumdmaps2

nc = 4;
nr = 4;
lgc = [0.5 0.5 0.5];
map_scale_fact = 0.3;
good_match = 0.0;
dxy = 7.5;
id = [16 15 12 11];

% Rotation fit, ipsilateral field
load VT1_rot_ipsi_field raxis maxscore lscore 
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_rotation_flow_map(raxis(1),raxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11 % Scale tmap to same lengths 
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
subplot(nr,nc,id(4))
hold off
for i = 1:11
    for j = 1:17
        if (lscore(i,j) > good_match)
            x1 = [az(1,j)-dxy, az(1,j)+dxy, az(1,j)+dxy, az(1,j)-dxy];
            y1 = [el(i,1)-dxy, el(i,1)-dxy, el(i,1)+dxy, el(i,1)+dxy];
            patch(x1, y1,lgc,'edgecolor','none')
            hold on
        end
    end
end
plot(raxis(1),raxis(2),'ok','markerfacecolor','k')
quiver(az,el,squeeze(tmap(1,:,:)),squeeze(tmap(2,:,:)),map_scale_fact,'color','r');
quiver(az,el,squeeze(cmap2(1,:,:)),squeeze(cmap2(2,:,:)),map_scale_fact,'k');
text(-135,100,['Rotation fit, ipsilateral field'])
text(-210,120,'B3','fontweight','bold','fontsize',14)

% Translation fit, ipsilateral field
load VT1_trans_ipsi_field taxis maxscore lscore 
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_translation_flow_map(taxis(1),taxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11 % Scale tmap to same lengths 
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
subplot(nr,nc,id(3))
hold off
for i = 1:11
    for j = 1:17
        if (lscore(i,j) > good_match)
            x1 = [az(1,j)-dxy, az(1,j)+dxy, az(1,j)+dxy, az(1,j)-dxy];
            y1 = [el(i,1)-dxy, el(i,1)-dxy, el(i,1)+dxy, el(i,1)+dxy];
            patch(x1, y1,lgc,'edgecolor','none')
            hold on
        end
    end
end
plot(taxis(1),taxis(2),'^k','markerfacecolor','k')
quiver(az,el,squeeze(tmap(1,:,:)),squeeze(tmap(2,:,:)),map_scale_fact,'color','r');
quiver(az,el,squeeze(cmap2(1,:,:)),squeeze(cmap2(2,:,:)),map_scale_fact,'k');
text(-135,100,['Translation fit, ipsilateral field'])
text(-210,120,'B4','fontweight','bold','fontsize',14)

% Translation fit, whole field
load VT1_trans_whole_field taxis maxscore lscore 
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_translation_flow_map(taxis(1),taxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11 % Scale tmap to same lengths 
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
subplot(nr,nc,id(1))
hold off
for i = 1:11
    for j = 1:17
        if (lscore(i,j) > good_match)
            x1 = [az(1,j)-dxy, az(1,j)+dxy, az(1,j)+dxy, az(1,j)-dxy];
            y1 = [el(i,1)-dxy, el(i,1)-dxy, el(i,1)+dxy, el(i,1)+dxy];
            patch(x1, y1,lgc,'edgecolor','none')
            hold on
        end
    end
end
plot(taxis(1),taxis(2),'^k','markerfacecolor','k')
quiver(az,el,squeeze(tmap(1,:,:)),squeeze(tmap(2,:,:)),map_scale_fact,'color','r');
quiver(az,el,squeeze(cmap2(1,:,:)),squeeze(cmap2(2,:,:)),map_scale_fact,'k');
text(-135,100,['Translation fit, Whole field'])
text(-210,120,'C4','fontweight','bold','fontsize',14)

% Rotation fit, whole field
load VT1_rot_whole_field raxis maxscore lscore 
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_rotation_flow_map(raxis(1),raxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11 % Scale tmap to same lengths 
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
subplot(nr,nc,id(2))
hold off
for i = 1:11
    for j = 1:17
        if (lscore(i,j) > good_match)
            x1 = [az(1,j)-dxy, az(1,j)+dxy, az(1,j)+dxy, az(1,j)-dxy];
            y1 = [el(i,1)-dxy, el(i,1)-dxy, el(i,1)+dxy, el(i,1)+dxy];
            patch(x1, y1,lgc,'edgecolor','none')
            hold on
        end
    end
end
raxis(1) = raxis(1) + 180;
plot(raxis(1),raxis(2),'ok','markerfacecolor','k')
quiver(az,el,squeeze(tmap(1,:,:)),squeeze(tmap(2,:,:)),map_scale_fact,'color','r');
quiver(az,el,squeeze(cmap2(1,:,:)),squeeze(cmap2(2,:,:)),map_scale_fact,'k');
text(-135,100,['Rotation fit, whole field'])
text(-210,120,'C3','fontweight','bold','fontsize',14)

for n = 1:4
    subplot(nr,nc,id(n))
    axis equal
    axis([-135 135 -90 90])
    set(gca,'tickdir','out')
    set(gca,'xtick',-135:15:135)
    set(gca,'ytick',-90:15:90)
    set(gca,'xticklabel',{'-135','','','-90','','','-45','','','0','','','45','','','90','','','135'})
    set(gca,'yticklabel',{'-90','','','-45','','','0','','','45','','','90'})
    xlabel('Azimuth (^o)')
    ylabel('Elevation (^o)')
end

