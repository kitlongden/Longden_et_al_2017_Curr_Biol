
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
    eval(['load ' datdir{nc} 'Map_data.mat x y u v']);
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

for n = 1:length(x)
    idx = find(ix(n,:) == x(n));
    idy = find(iy(:,n) == y(n));
    imu(idy,idx) = mu_RF(n);
    imv(idy,idx) = mv_RF(n);
end

cmap = zeros(2,11,17);
cmap(1,:,:) = imu;
cmap(2,:,:) = imv;

az = ix;
el = iy;

cmap2 = cmap;
for i = 1:11
    for j = 1:17 % Set amplitude of cell map to unity
        cmap2(:,i,j) = cmap2(:,i,j)/sqrt(cmap(1,i,j)^2 + cmap(2,i,j)^2);
    end
end

camp = zeros(11,17);
for i = 1:11
    for j = 1:17
        camp(i,j) = sqrt(cmap(1,i,j)^2 + cmap(2,i,j)^2);
    end
end

% Summed dot product of cell map
dmaps2 = squeeze(dot(cmap,cmap)); % data is 11 x 17
azsumdmaps2 = sum(dmaps2,2); % data i2 11
sumdmaps2 = sum(azsumdmaps2.*cosd(el(:,1)));

save Figure_1_mean_data.mat az el cmap cmap2 camp sumdmaps2
    
    