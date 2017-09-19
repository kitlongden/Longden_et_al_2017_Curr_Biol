
addpath mapfit_fns/

%% Load mean data, calculated in part_2
% cmap = Mean data in figure 1 (cell map)
% cmap2 = cmap normalised to vectors of unit length
% camp = array of vector lengths in cmap
% sumdmaps = sum of dot products of cmap with itself
load Figure_1_mean_data.mat az el cmap cmap2 camp sumdmaps2

%% Translation fit, whole field

% Search parameters
dang = 1;
az_ang = -180:dang:180;
el_ang = -90:dang:90;
N_az = length(az_ang);
N_el = length(el_ang);

% Fits for direction only


% Set up search
maxscore = 0;
taxis = [-1000,-1000];
score = zeros(N_az,N_el);

for na = 1:N_az
    if (mod(na,10) == 0)
        disp(num2str(na))
    end
    for ne = 1:N_el
        tmap = zeros(2,11,17);
        [az,el,u,v] = Generate_translation_flow_map(az_ang(na),el_ang(ne),.1);
        tmap(1,:,:) = u;
        tmap(2,:,:) = v;
        
        % Scale tmap to unit lengths
        for i = 1:11
            for j = 1:17
                tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
            end
        end
        
        % Integrate dot products, was cos elevation sampling correction
        dmaps = squeeze(dot(tmap,cmap2)); % data is 11 x 17
        azsumdmaps = sum(dmaps,2); % data i2 11
        sumdmaps = sum(azsumdmaps.*cosd(el(:,1)));
        
        score(na,ne) = sumdmaps/(11*17);
        if (score(na,ne) > maxscore)
            taxis = [az_ang(na), el_ang(ne)];
            maxscore = score(na,ne);
        end
    end
end

% Local scores for best fit
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_translation_flow_map(taxis(1),taxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
lscore = zeros(11,17);
for i = 1:11
    for j= 1:17
        lscore(i,j) = dot(tmap(:,i,j),cmap2(:,i,j))/dot(cmap2(:,i,j),cmap2(:,i,j));
    end
end

save VT1_trans_whole_field taxis maxscore score az_ang el_ang lscore

%% Rotation fit, whole field

% Set up search
maxscore = 0;
raxis = [-1000,-1000];
score = zeros(N_az,N_el);

for na = 1:N_az
    if (mod(na,10) == 0)
        disp(num2str(na))
    end
    for ne = 1:N_el
        rmap = zeros(2,11,17);
        [az,el,u,v] = Generate_rotation_flow_map(az_ang(na),el_ang(ne),.1);
        rmap(1,:,:) = u;
        rmap(2,:,:) = v;
        
        % Scale tmap to unit lengths
        for i = 1:11
            for j = 1:17
                rmap(:,i,j) = rmap(:,i,j)/sqrt(rmap(1,i,j)^2 + rmap(2,i,j)^2);
            end
        end
        
        % Integrate dot products, was cos elevation sampling correction
        dmaps = squeeze(dot(rmap,cmap2)); % data is 11 x 17
        azsumdmaps = sum(dmaps,2); % data i2 11
        sumdmaps = sum(azsumdmaps.*cosd(el(:,1)));
        
        score(na,ne) = sumdmaps/(11*17);
        if (score(na,ne) > maxscore)
            raxis = [az_ang(na), el_ang(ne)];
            maxscore = score(na,ne);
        end
    end
end

% Local scores for best fit
rmap = zeros(2,11,17);
[az,el,u,v] = Generate_rotation_flow_map(raxis(1),raxis(2),.1);
rmap(1,:,:) = u;
rmap(2,:,:) = v;
for i = 1:11
    for j = 1:17
        rmap(:,i,j) = rmap(:,i,j)/sqrt(rmap(1,i,j)^2 + rmap(2,i,j)^2);
    end
end
lscore = zeros(11,17);
for i = 1:11
    for j= 1:17
        lscore(i,j) = dot(rmap(:,i,j),cmap2(:,i,j))/dot(cmap2(:,i,j),cmap2(:,i,j));
    end
end

save VT1_rot_whole_field raxis maxscore score az_ang el_ang lscore

%% Translation fit, ipsilateral field

% Set up search
maxscore = 0;
taxis = [-1000,-1000];
score = zeros(N_az,N_el);

% Azimuth window for defining ipsilateral responses
az_start = 9;
az_end = 17;

for na = 1:N_az
    if (mod(na,10) == 0)
        disp(num2str(na))
    end
    for ne = 1:N_el
        tmap = zeros(2,11,17);
        [az,el,u,v] = Generate_translation_flow_map(az_ang(na),el_ang(ne),.1);
        tmap(1,:,:) = u;
        tmap(2,:,:) = v;
        
        % Scale tmap to unit lengths
        for i = 1:11
            for j = 1:17
                tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
            end
        end
        
        % Integrate dot products, was cos elevation sampling correction
        dmaps = squeeze(dot(tmap(:,:,az_start:az_end),cmap2(:,:,az_start:az_end))); % data is 11 x 17
        azsumdmaps = sum(dmaps,2); % data i2 11
        sumdmaps = sum(azsumdmaps.*cosd(el(:,1)));
        
        score(na,ne) = sumdmaps/(11*(length(az_start:az_end)));
        if (score(na,ne) > maxscore)
            taxis = [az_ang(na), el_ang(ne)];
            maxscore = score(na,ne);
        end
    end
end

% Local scores for best fit
tmap = zeros(2,11,17);
[az,el,u,v] = Generate_translation_flow_map(taxis(1),taxis(2),.1);
tmap(1,:,:) = u;
tmap(2,:,:) = v;
for i = 1:11
    for j = 1:17
        tmap(:,i,j) = tmap(:,i,j)/sqrt(tmap(1,i,j)^2 + tmap(2,i,j)^2);
    end
end
lscore = zeros(11,17);
for i = 1:11
    for j= 1:17
        lscore(i,j) = dot(tmap(:,i,j),cmap2(:,i,j))/dot(cmap2(:,i,j),cmap2(:,i,j));
    end
end

save VT1_trans_ipsi_field taxis maxscore score az_ang el_ang lscore

%% Rotation fit, ipsilateral field

% Azimuth window for defining ipsilateral responses
az_start = 9;
az_end = 17;

% Set up search
maxscore = 0;
raxis = [-1000,-1000];
score = zeros(N_az,N_el);

for na = 1:N_az
    if (mod(na,10) == 0)
        disp(num2str(na))
    end
    for ne = 1:N_el
        rmap = zeros(2,11,17);
        [az,el,u,v] = Generate_rotation_flow_map(az_ang(na),el_ang(ne),.1);
        rmap(1,:,:) = u;
        rmap(2,:,:) = v;
        
        % Scale tmap to unit lengths
        for i = 1:11
            for j = 1:17
                rmap(:,i,j) = rmap(:,i,j)/sqrt(rmap(1,i,j)^2 + rmap(2,i,j)^2);
            end
        end
        
        % Integrate dot products, was cos elevation sampling correction
        dmaps = squeeze(dot(rmap(:,:,az_start:az_end),cmap2(:,:,az_start:az_end))); % data is 11 x 17
        azsumdmaps = sum(dmaps,2); % data i2 11
        sumdmaps = sum(azsumdmaps.*cosd(el(:,1)));
        
        score(na,ne) = sumdmaps/(11*length(az_start:az_end));
        if (score(na,ne) > maxscore)
            raxis = [az_ang(na), el_ang(ne)];
            maxscore = score(na,ne);
        end
    end
end
% Local scores for best fit
rmap = zeros(2,11,17);
[az,el,u,v] = Generate_rotation_flow_map(raxis(1),raxis(2),.1);
rmap(1,:,:) = u;
rmap(2,:,:) = v;
for i = 1:11
    for j = 1:17
        rmap(:,i,j) = rmap(:,i,j)/sqrt(rmap(1,i,j)^2 + rmap(2,i,j)^2);
    end
end
lscore = zeros(11,17);
for i = 1:11
    for j= 1:17
        lscore(i,j) = dot(rmap(:,i,j),cmap2(:,i,j))/dot(cmap2(:,i,j),cmap2(:,i,j));
    end
end

save VT1_rot_ipsi_field raxis maxscore score az_ang el_ang lscore
