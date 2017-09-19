

export_flag = 0;
nx = 4;
ny = 8;


%%
bucol = zeros(3,3);
bucol(1,:) = [0 0 1];
bucol(2,:) = [0.5 0 0.5];
bucol(3,:) = [1 0 0];
lgc = [0.5 0.5 0.5];
lgc2 = 0.75*[1 1 1];
mks = 5;


h_fig4 = figure(4);
set(h_fig4,'color','w','Position',[100 200 800 600]) % 800 wide for 17.6 cm (2 columns, J Neuro)

%%

N_looms = 8;
t = (1:80000)/20000;
load r_burst_data.mat
afbarray_1(:,end-400:end) = NaN;
afbarray_2(:,end-400:end) = NaN;
afbarray_3(:,end-400:end) = NaN;
afbarray_1(:,1:1600) = NaN;
afbarray_2(:,1:1600) = NaN;
afbarray_3(:,1:1600) = NaN;

for nl = 1:N_looms/2
    subplot(ny,nx,[nl+nx nl+2*nx])
    hold off
    plot(t,afbarray_3(nl,:),'color',bucol(3,:));
    hold on
    plot(t,afbarray_1(nl,:),'color',bucol(1,:));

    box off
    set(gca,'tickdir','out')
    set(gca,'xtick',0:0.5:4,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'})
    axis([0 4 0 25])
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    if (nl == 1)
        legend('3+ spike bursts', 'Single spikes')
        legend boxoff
        text(-0.75,25,'B','fontweight','bold','fontsize',14)
    end
end
for nl = 1+N_looms/2:N_looms
    subplot(ny,nx,[nl+3*nx-4 nl+4*nx-4])
    hold off
    plot(t,afbarray_3(nl,:),'color',bucol(3,:));
    hold on
    plot(t,afbarray_1(nl,:),'color',bucol(1,:));

    box off
    set(gca,'tickdir','out')
    set(gca,'xtick',0:0.5:4,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'})
    axis([0 4 0 25])
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    if (nl == 5)
        text(-0.75,25,'C','fontweight','bold','fontsize',14)
    end
end

%% Plot stimuli

screen_dist = 200; % mm
screen_height = 263; % mm
screen_width = 359; % mm
loom_object_radius = 50; % mm
loom_duration = 3; % seconds
loom_radius_max = 41; % degrees. 
% This is because: tand(theta) = (screen_width/2)/screen_distance |-> theta = 41
% So the screen subtends 82 degrees and the width of the loom object
% subtends that at the end.
loom_radius_max_distance = loom_object_radius / tand(loom_radius_max);
loom_speeds = [30 150 500 5000]; % mm/s
r_over_v = (loom_object_radius./loom_speeds) *1000; % ms
N_looms = length(loom_speeds);

dt = 1/200;
t2 = dt:dt:4;
N_t = length(t2);
stim = zeros(N_looms, N_t);
dist = zeros(N_looms, N_t);
stim2 = zeros(N_looms, N_t); % for calculating trajectory of stimuli 
dist2 = zeros(N_looms, N_t); % if they had continued.
t_collision = zeros(N_looms,1);
for nl = 1:N_looms
    dx = loom_speeds(nl)/200;
    dist(nl,101:700) =  loom_radius_max_distance + dx*599 : -dx : loom_radius_max_distance;
    dist(nl,701:800) = loom_radius_max_distance;
    stim(nl,101:800) = 2*atand(50./dist(nl,101:800));
    
    dist2(nl,:) = dist(nl,:);
    dist2(nl,701:800) = loom_radius_max_distance - dx: -dx : loom_radius_max_distance - dx*100;
    id = 100 + find(dist2(nl,101:800) <= 0, 1,'first');
    if (isempty(id))
        t_collision(nl) = NaN;
    else
        t_collision(nl) = id/200;
    end
    % dist2(nl,id(nl):800) = 0;
    stim2(nl,101:800) = 2*atand(50./dist2(nl,101:800));
end

for nl = 1:N_looms
    subplot(ny,nx,nl)
    hold off
    plot(t2,stim(nl,:),'color',lgc)
    hold on
    line([t_collision(1) t_collision(1)], [0 90], 'color', lgc)
    line([t_collision(2) t_collision(2)], [0 90], 'color', lgc2)
    line([t_collision(4) t_collision(4)], [0 90], 'color', 'k')
    axis([0 4 0 90])
    box off
    set(gca,'tickdir','out')
    set(gca,'xtick',0:0.5:4,'xticklabel',{'0' '' '' '' '2' '' '3' '3.5' '4'})
    set(gca,'ytick',0:30:90,'yticklabel',{'0' '' '' '90'});
    xlabel('Time (s)')
    ylabel('Angle')
    if (nl < 5)
        text(0.5,80,['r / v = ' num2str(r_over_v(nl)) ' (ms)'])
    else
        text(0.5,80,['r / v = ' num2str(r_over_v(nl-4)) ' (ms)'])
    end
    %     legend('3', '15', '500','location','northwest')
    %     legend boxoff
    if (nl == 1)
        text(-0.75,90,'A','fontweight','bold','fontsize',14)
    end
end

%%
nx = 6;

N_cells = 11;
% Stage3 is t=3-3.5
mnrstage3_1 = mean(r_stage3_1,3);
mnrstage3_3 = mean(r_stage3_3,3);
mn_mrstage3_1 = mean(mnrstage3_1,1);
se_mrstage3_1 = std(mnrstage3_1,1)/sqrt(N_cells);
mn_mrstage3_3 = mean(mnrstage3_3,1);
se_mrstage3_3 = std(mnrstage3_3,1)/sqrt(N_cells);

% Stage4 is t=3.5-4
mnrstage4_1 = mean(r_stage4_1,3);
mnrstage4_3 = mean(r_stage4_3,3);
mn_mrstage4_1 = mean(mnrstage4_1,1);
se_mrstage4_1 = std(mnrstage4_1,1)/sqrt(N_cells);
mn_mrstage4_3 = mean(mnrstage4_3,1);
se_mrstage4_3 = std(mnrstage4_3,1)/sqrt(N_cells);

subplot(ny,nx,[31 32 43 44])
hold off
semilogx(r_over_v,mn_mrstage3_3(1:4),'-^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
hold on
semilogx(r_over_v,mn_mrstage3_1(1:4),'-o','color',bucol(1,:),'markerfacecolor',bucol(1,:))
semilogx(r_over_v,mn_mrstage4_3(1:4),'--^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
semilogx(r_over_v,mn_mrstage4_1(1:4),'--o','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(r_over_v,mn_mrstage3_3(1:4),se_mrstage3_3(1:4),'-^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
errorbar(r_over_v,mn_mrstage3_1(1:4),se_mrstage3_3(1:4),'-o','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(r_over_v,mn_mrstage3_3(5:8),se_mrstage3_3(5:8),'--^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
errorbar(r_over_v,mn_mrstage3_1(5:8),se_mrstage3_3(5:8),'--o','color',bucol(1,:),'markerfacecolor',bucol(1,:))

box off
axis([5 3000 0 30])
set(gca,'tickdir','out')
xlabel('r/v (ms)')
ylabel('Spike rate (Hz)')
legend('3+ spike bursts, dark loom', 'Single spikes, dark loom', '3+ spike bursts, bullseye', 'Single spikes, bullseye')
legend boxoff
text(10,25,{'Mean spike' 'rate 3-3.5 s'})
set(gca,'xtick',[5:10 20:10:100 200:100:1000 2000 3000],'xticklabel',{'' '' '' '' '' '10' '' '' '' '' '' '' '' '' '100' '' '' '' '' '' '' '' '' '1000' '' ''})
text(4,30,'E','fontweight','bold','fontsize',14)

subplot(ny,nx,[33 34 45 46])
hold off
semilogx(r_over_v,mn_mrstage4_3(1:4),'-^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
hold on
semilogx(r_over_v,mn_mrstage4_1(1:4),'-o','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(r_over_v,mn_mrstage4_3(1:4),se_mrstage4_3(1:4),'-^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
errorbar(r_over_v,mn_mrstage4_1(1:4),se_mrstage4_3(1:4),'-o','color',bucol(1,:),'markerfacecolor',bucol(1,:))
errorbar(r_over_v,mn_mrstage4_3(5:8),se_mrstage4_3(5:8),'--^','color',bucol(3,:),'markerfacecolor',bucol(3,:))
errorbar(r_over_v,mn_mrstage4_1(5:8),se_mrstage4_3(5:8),'--o','color',bucol(1,:),'markerfacecolor',bucol(1,:))

box off
axis([5 3000 0 30])
set(gca,'tickdir','out')
xlabel('r/v (ms)')
ylabel('Spike rate (Hz)')
text(10,25,{'Mean spike' 'rate 3.5-4 s'})
set(gca,'xtick',[5:10 20:10:100 200:100:1000 2000 3000],'xticklabel',{'' '' '' '' '' '10' '' '' '' '' '' '' '' '' '100' '' '' '' '' '' '' '' '' '1000' '' ''})











