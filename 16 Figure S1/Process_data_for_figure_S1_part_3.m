function Process_data_for_figure_S1_part_3()

addpath mapfit_fns/

datadir = '../Data/''Local motion receptive fields''/VT1/Wide'' field maps''/';
fname = 'Map_data.mat';
fdir{1} = '16'' Apr 2008''/';
fdir{2} = '21'' Mar 2007''/';
fdir{3} = '29'' Apr 2008''/';
fdir{4} = '08'' Jan 2014-1''/';
fdir{5} = '09'' Jan 2014-1''/';
fdir{6} = '09'' Jan 2014-2''/';
fdir{7} = '10'' Jan 2014-1''/';
fdir{8} = '14'' Jan 2014-1''/';
fdir{9} = '21'' Jan 2014-3''/';
fdir{10} = '03'' Dec 2002''/';
N_cells = length(fdir);

N_con = 4;
axis = zeros(N_con,N_cells,2);
score = zeros(N_con,N_cells,1);

% Set up cell map data
for nc = 1:N_cells
    setup_cmaps(nc,datadir,fdir,fname);
end

% Fits
for nc = 1:N_cells
    disp(['Cell ...' num2str(nc)])
    
    disp('...translation...')
    trans_rot_flag = 1;  % trans_rot_flag = 1 for translation; trans_rot_flag ~= 1 for rotation
    flow_fit(nc,trans_rot_flag);
    
    disp('...rotatation...')
    trans_rot_flag = 0;
    flow_fit(nc,trans_rot_flag);
end

end

function setup_cmaps(nc,datadir,fdir,fname)

[az,el] = meshgrid(-120:15:120,75:-15:-75);
el(11,:) = -70;

% load data
if (nc <= 9)
    eval(['load ' datadir fdir{nc} fname]);
else
    eval(['load ' datadir fdir{nc} 'dot_mapping1_03-Dec-2002-SJH/Data_LPDS_xyuv.mat']);
end
% cellmap
cmap = zeros(2,11,17);
if (length(u) == 66)
    c = interpmap(u,v);
    cmap(1,:,:) = c{1};
    cmap(2,:,:) = c{2};
else
    cmap(1,:,:) = u;
    cmap(2,:,:) = v;
end
% Set amplitude of cell map to unity
cmap2 = cmap;
for i = 1:11
    for j = 1:17
        if (sqrt(cmap(1,i,j)^2 + cmap(2,i,j)^2) > 0)
            cmap2(:,i,j) = cmap2(:,i,j)/sqrt(cmap(1,i,j)^2 + cmap(2,i,j)^2);
        else
            cmap2(:,i,j) = 0;
        end
    end
end
% Amplitude of cell map
camp = zeros(11,17);
for i = 1:11
    for j = 1:17
        camp(i,j) = sqrt(cmap(1,i,j)^2 + cmap(2,i,j)^2);
    end
end

eval(['save Cell_maps/cell_' num2str(nc) '_cmaps.mat az el cmap cmap2 camp'])

end

function flow_fit(nc,trans_rot_flag)

% Search parameters
dang = 1;
az_ang = -180:dang:180;
el_ang = -90:dang:90;
N_az = length(az_ang);
N_el = length(el_ang);
maxscore = 0;
fscore = zeros(N_az,N_el);
axes = [-1000,-1000];

% Azimuth window for defining ipsilateral responses
az_start = 1;
az_end = 17;

% Load data
eval(['load Cell_maps/cell_' num2str(nc) '_cmaps.mat el cmap cmap2 camp'])

% Summed dot product of cell map
dmaps2 = squeeze(dot(cmap(:,:,az_start:az_end),cmap(:,:,az_start:az_end))); % data is 11 x 17
azsumdmaps2 = sum(dmaps2,2); % data i2 11
sumdmaps2 = sum(azsumdmaps2.*cosd(el(:,1)));

% Summed dot product of normalized cell map
dmaps3 = squeeze(dot(cmap2(:,:,az_start:az_end),cmap2(:,:,az_start:az_end))); % data is 11 x 17
azsumdmaps3 = sum(dmaps3,2); % data i2 11
sumdmaps3 = sum(azsumdmaps3.*cosd(el(:,1)));

display('...searching angle...')
for na = 1:N_az
    if (mod(na,60) == 0)
        disp(['...' num2str(na)])
    end
    for ne = 1:N_el
        map = zeros(2,11,17);
        if (trans_rot_flag == 1)
            [az,el,u,v] = Generate_translation_flow_map(az_ang(na),el_ang(ne),.1);
        else
            [az,el,u,v] = Generate_rotation_flow_map(az_ang(na),el_ang(ne),.1);
        end
        map(1,:,:) = u;
        map(2,:,:) = v;
        
        % Scale tmap to unit lengths
        for i = 1:11
            for j = 1:17
                map(:,i,j) = map(:,i,j)/sqrt(map(1,i,j)^2 + map(2,i,j)^2);
            end
        end
        
        % Integrate dot products, was cos elevation sampling correction
        dmaps = squeeze(dot(map(:,:,az_start:az_end),cmap2(:,:,az_start:az_end))); % data was 11 x 17
        azsumdmaps = sum(dmaps,2); % data i2 11
        sumdmaps = sum(azsumdmaps.*cosd(el(:,1)));
        
        % Calculate score
        fscore(na,ne) = sumdmaps/sumdmaps3;
        
        % Identify best fit
        if (fscore(na,ne) > maxscore)
            axes = [az_ang(na), el_ang(ne)];
            maxscore = fscore(na,ne);
        end
    end
end

if (trans_rot_flag == 1)
    eval(['save Whole_field_cell_fits/Translation/cell_' num2str(nc) '.mat axes maxscore'])
else
    eval(['save Whole_field_cell_fits/Rotation/cell_' num2str(nc) '.mat axes maxscore'])
end

end

