

%%
nx = 4;
ny = 5;

%%
bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];
lgc = [0.5 0.5 0.5];
lgc2 = 0.75*[1 1 1];
mks = 5;

%%
N_step_heights = 8;
sh = [1 6 11 16 21 26 36 46]-26; % cm 26 is ground height
sh = sh/20;
sh2 = sh;
N_cells = 10;

h_fig5 = figure(5);
set(h_fig5,'color','w','Position',[100 200 1000 800]) % 800 wide for 17.6 cm (2 columns, J Neuro)

% Calculate spatial frequency tuning sine
load ../'07 Multiple depth planes'/r_burst_data_sine.mat
sin_mrstage3_1 = mean(r_stage3_1,3);
sin_mrstage3_2 = mean(r_stage3_2,3);
sin_mrstage3_3  = mean(r_stage3_3,3);
sin_mrstage3_12 = mean(r_stage3_1+r_stage3_2,3);
sin_mrstage3_123 = mean(r_stage3_1+r_stage3_2+r_stage3_3,3);
sin_mrstage4_1 = mean(r_stage4_1,3);
sin_mrstage4_2 = mean(r_stage4_2,3);
sin_mrstage4_3  = mean(r_stage4_3,3);
sin_mrstage4_12 = mean(r_stage4_1+r_stage4_2,3);
sin_mrstage4_123 = mean(r_stage4_1+r_stage4_2+r_stage4_3,3);

% Calculate spatial frequency tuning grey
load ../'07 Multiple depth planes'/r_burst_data_grey.mat
grey_mrstage3_1 = mean(r_stage3_1,3);
grey_mrstage3_2 = mean(r_stage3_2,3);
grey_mrstage3_3  = mean(r_stage3_3,3);
grey_mrstage3_12 = mean(r_stage3_1+r_stage3_2,3);
grey_mrstage3_123 = mean(r_stage3_1+r_stage3_2+r_stage3_3,3);
grey_mrstage4_1 = mean(r_stage4_1,3);
grey_mrstage4_2 = mean(r_stage4_2,3);
grey_mrstage4_3  = mean(r_stage4_3,3);
grey_mrstage4_12 = mean(r_stage4_1+r_stage4_2,3);
grey_mrstage4_123 = mean(r_stage4_1+r_stage4_2+r_stage4_3,3);

% Plot sin
subplot(ny,nx,[3*nx+1 4*nx+1])
hold off
plot(sh,mean(sin_mrstage3_3), '^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1,'markersize',mks)
hold on
plot(sh,mean(sin_mrstage3_1),'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(sin_mrstage3_3),std(sin_mrstage3_3)/sqrt(N_cells),'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(sin_mrstage3_1),std(sin_mrstage3_1)/sqrt(N_cells),'s-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
text(-1.25,49,'Sine background','fontsize',12,'color','k','fontweight','bold')
text(-1.25,45,'t = 2-3 s','fontsize',12,'color','k')

box off
set(gca,'tickdir','out')
set(gca,'xtick',sh2,'xticklabel',{'' '-1' '' '' '' '0' '' '1'})
set(gca,'ytick',0:10:50,'yticklabel',{'0' '10' '20' '30' '40' ''})
axis([-1.45 1.2 0 50])
ylabel('Spike rate (Hz)')
xlabel('Stripe height (\lambda)')
text(-1.9,50,'F','fontweight','bold','fontsize',14)

legend('3+ spike bursts', 'single spikes','location','southeast')
legend boxoff

subplot(ny,nx,[3*nx+2 4*nx+2])
hold off
plot(sh,mean(sin_mrstage4_3), '^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
hold on
plot(sh,mean(sin_mrstage4_1),'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(sin_mrstage4_3),std(sin_mrstage4_3)/sqrt(N_cells),'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(sin_mrstage4_1),std(sin_mrstage4_1)/sqrt(N_cells),'s-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
text(-1.25,49,'Sine background','fontsize',12,'color','k','fontweight','bold')
text(-1.25,45,'t = 3-3.5 s','fontsize',12,'color','k')

box off
set(gca,'tickdir','out')
set(gca,'xtick',sh2,'xticklabel',{'' '-1' '' '' '' '0' '' '1'})
set(gca,'ytick',0:10:50,'yticklabel',{'0' '10' '20' '30' '40' ''})
axis([-1.45 1.2 0 50])
ylabel('Spike rate (Hz)')
xlabel('Stripe height (\lambda)')
text(-1.9,50,'G','fontweight','bold','fontsize',14)

% Plot grey
subplot(ny,nx,[3*nx+3 4*nx+3])
hold off
plot(sh,mean(grey_mrstage3_3), '^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
hold on
plot(sh,mean(grey_mrstage3_1),'s-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(grey_mrstage3_3),std(grey_mrstage3_3)/sqrt(N_cells),'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(grey_mrstage3_1),std(grey_mrstage3_1)/sqrt(N_cells),'s-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
text(-1.25,49,'Gray background','fontsize',12,'color','k','fontweight','bold')
text(-1.25,45,'t = 2-3 s','fontsize',12,'color','k')

box off
set(gca,'tickdir','out')
set(gca,'xtick',sh2,'xticklabel',{'' '-1' '' '' '' '0' '' '1'})
set(gca,'ytick',0:10:50,'yticklabel',{'0' '10' '20' '30' '40' ''})
axis([-1.45 1.2 0 50])
ylabel('Spike rate (Hz)')
xlabel('Stripe height (\lambda)')
text(-1.9,50,'H','fontweight','bold','fontsize',14)

subplot(ny,nx,[3*nx+4 4*nx+4])
hold off
plot(sh,mean(grey_mrstage4_3), '^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
hold on
plot(sh,mean(grey_mrstage4_1),'o-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(grey_mrstage4_3),std(grey_mrstage4_3)/sqrt(N_cells),'^-','color',bucol(3,:),'markerfacecolor',bucol(3,:),'linewidth',1.0,'markersize',mks)
errorbar(sh,mean(grey_mrstage4_1),std(grey_mrstage4_1)/sqrt(N_cells),'s-','color',bucol(1,:),'markerfacecolor',bucol(1,:),'linewidth',1.0,'markersize',mks)
text(-1.25,49,'Gray background','fontsize',12,'color','k','fontweight','bold')
text(-1.25,45,'t = 3-3.5 s','fontsize',12,'color','k')

box off
set(gca,'tickdir','out')
set(gca,'xtick',sh2,'xticklabel',{'' '-1' '' '' '' '0' '' '1'})
set(gca,'ytick',0:10:50,'yticklabel',{'0' '10' '20' '30' '40' ''})
axis([-1.45 1.2 0 50])
ylabel('Spike rate (Hz)')
xlabel('Stripe height (\lambda)')
text(-1.9,50,'J','fontweight','bold','fontsize',14)



%%
nx = 16;
ny = 5;
did = 4;
id = zeros(8,2);
id(1,:) = [1 3];
id(2,:) = [1 3]; %
id(3,:) = [1 3];
id(4,:) = [4 6];
id(5,:) = [7 9];
id(6,:) = [7 9];
id(8,:) = [10 12];

N_step_heights = 8;
t = (1:80000)/20000;

% Remove artefacts of filtering
load ../'07 Multiple depth planes'/r_burst_data_sine.mat
afbarray_1(:,1:1400) = NaN;
afbarray_2(:,1:1400) = NaN;
afbarray_3(:,1:1400) = NaN;
afbarray_1(:,end-400:end) = NaN;
afbarray_2(:,end-400:end) = NaN;
afbarray_3(:,end-400:end) = NaN;

for nsh = [3 4 6 8]
    subplot(ny,nx,nx+did+id(nsh,:))
    hold off
    plot(t,afbarray_3(nsh,:),'color',bucol(3,:));
    hold on
    plot(t,afbarray_1(nsh,:),'color',bucol(1,:));
    
    box off
    set(gca,'tickdir','out')
    set(gca,'xtick',0:0.5:4,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'})
    axis([0 4 0 75])
    xlabel('Time (s)')
    if (nsh == 3)
        ylabel('Spike rate (Hz)')
        set(gca,'ytick',0:10:65,'yticklabel',{'0' '' '20' '' '40' '' '60'})
        legend('3+ spike bursts', 'single spikes')
        legend boxoff
    else
        set(gca,'ytick',0:10:65,'yticklabel',{'' '' '' '' '' '' ''})
    end
end

% Remove artefacts of filtering
load ../'07 Multiple depth planes'/r_burst_data_grey.mat
afbarray_1(:,1:1400) = NaN;
afbarray_2(:,1:1400) = NaN;
afbarray_3(:,1:1400) = NaN;
afbarray_1(:,end-400:end) = NaN;
afbarray_2(:,end-400:end) = NaN;
afbarray_3(:,end-400:end) = NaN;

for nsh = [3 4 6 8]
    subplot(ny,nx,2*nx+did+id(nsh,:))
    hold off
    plot(t,afbarray_3(nsh,:),'color',bucol(3,:));
    hold on
    plot(t,afbarray_1(nsh,:),'color',bucol(1,:));
    
    box off
    set(gca,'tickdir','out')
    set(gca,'xtick',0:0.5:4,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'})
    axis([0 4 0 75])
    xlabel('Time (s)')
    if (nsh == 3)
        ylabel('Spike rate (Hz)')
        set(gca,'ytick',0:10:65,'yticklabel',{'0' '' '20' '' '40' '' '60'})
    else
        set(gca,'ytick',0:10:65,'yticklabel',{'' '' '' '' '' '' ''})
    end
end

subplot(ny,nx,nx+did+id(3,:))
fs = 10;
text(-0.75,75,'C','fontweight','bold','fontsize',14)
text(0.1,55,'Stripe height = -3\lambda/4','fontsize',fs)
text(0.1,70,'Sine background','fontsize',fs,'fontweight','bold')
subplot(ny,nx,nx+did+id(4,:))
text(0.1,55,'Stripe height = -\lambda/2','fontsize',fs)
subplot(ny,nx,nx+did+id(6,:))
text(0.1,55,'Stripe height = 0','fontsize',fs)
subplot(ny,nx,nx+did+id(8,:))
text(0.1,55,'Stripe height = \lambda','fontsize',fs)


subplot(ny,nx,2*nx+did+id(3,:))
text(-0.95,75,'D','fontweight','bold','fontsize',14)
text(0.1,55,'Stripe height = -3\lambda/4','fontsize', fs)
text(0.1,70,'Gray background','fontsize',fs,'fontweight','bold')
subplot(ny,nx,2*nx+did+id(4,:))
text(0.1,55,'Stripe height = -\lambda/2','fontsize',fs)
subplot(ny,nx,2*nx+did+id(6,:))
text(0.1,55,'Stripe height = 0' ,'fontsize',fs)
subplot(ny,nx,2*nx+did+id(8,:))
text(0.1,55,'Stripe height = \lambda','fontsize',fs)

%% Plot stimuli
N_frames = 800;
load ../'07 Multiple depth planes'/spacetimeplots_sine.mat ground
ground_sine = ground;
load ../'07 Multiple depth planes'/spacetimeplots_grey.mat;
ground_grey = ground;

for nsh = [3 4 6 8]
    
    subplot(ny,nx,4+id(nsh,:))
    hold off
    mesh(ground_sine{nsh})
    axis([0 N_frames 0 480 0 2])
    view(2)
    colormap('bone')
    set(gca,'xtick',0:100:N_frames,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'});
    % axis off
    xlabel('Time (s)')
    if (nsh == 3)
        ylabel('Elevation')
        set(gca,'ytick',[0 480],'yticklabel',{'-52^o' '0'})
    else
        set(gca,'ytick',[0 480],'yticklabel',{'' ''})
    end
    set(gca,'tickdir','out')
    
end
subplot(ny,nx,4+id(3,:))
text(-150, 500, 'B','fontweight','bold','fontsize',14)
text(0,525,'Spacetime plots of visual stimulus')
