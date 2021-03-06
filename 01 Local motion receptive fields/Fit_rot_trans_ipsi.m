function response = Fit_rot_trans_ipsi(rotation,normscheme,directionscheme,azstart,elstart,u,v,Ni,Nj,figureflag)
% If <rotation=1> then the code fits a rotation, if not then it fits a translation

% Azimuth and elevation for interpolated maps
[az,el] = meshgrid(-120:15:120,75:-15:-75);
el(11,:) = -70;

% Get the x and y components of the LPDS
cellmap{1} = zeros(11,17);
cellmap{2} = zeros(11,17);
if (length(u) == 66)
    cellmap = interpmap(u,v);
else
    cellmap{1} = u;
    cellmap{2} = v;
end

% Normalise maps 
% ...so the mean LMS value is 1, (normscheme = 1), 
% ...or by the max LMS value (normscheme = 2),
% ...or for equal magnitude arrows (normscheme = 3).
anorm = 0;
atemp = zeros(11,17);
atemp(:,:) = sqrt(cellmap{1}.^2 + cellmap{2}.^2);
if (normscheme == 1)
    anorm = mean(mean(atemp));
    cellmap{1} = cellmap{1}/anorm;
    cellmap{2} = cellmap{2}/anorm;
elseif (normscheme == 2)
    anorm = max(max(atemp));
    cellmap{1} = cellmap{1}/anorm;
    cellmap{2} = cellmap{2}/anorm;
elseif (normscheme == 3)
    for ne = 1:11
        for na = 1:17
            if (atemp(ne,na) > 0)
                cellmap{1}(ne,na) = cellmap{1}(ne,na)/atemp(ne,na);
                cellmap{2}(ne,na) = cellmap{2}(ne,na)/atemp(ne,na);
            else
                cellmap{1}(ne,na) = 0;
                cellmap{2}(ne,na) = 0;
            end
        end
    end
end

cmap = zeros(2,11,17);
cmap(1,:,:) = cellmap{1};
cmap(2,:,:) = cellmap{2};


% Set up match parameters
maxscore = 0;
maxaxis = [-360,-360];
score = zeros(Ni,Nj);

% loop through azimuth and elevation angles for the axis of rotation,
% and elevation factor (elevationfactor = 1)
for i = 1:Ni
    % disp(Ni-i);
    for j = 1:Nj
        azi = azstart + 1*(i-1);
        elj = elstart + 1*(j-1);

        % Get a rotation or translation flow map
        % N.b. rotaz and rotel are the x and y components of the flow map
        if (rotation == 1)
            flowmap = Generate_rotation_flow_map(azi,elj,0.1);
        else
            flowmap = Generate_translation_flow_map(azi,elj,0.1);
        end
        
        % Scale the flow map so that the magnitudes match
        for i2 = 1:11
            for j2 = 1:17
                length_cmap = sqrt(cmap(1,i2,j2).^2 + cmap(2,i2,j2).^2);
                length_fmap = sqrt(flowmap(1,i2,j2).^2 + flowmap(2,i2,j2).^2);
                flowmap(:,i2,j2) = flowmap(:,i2,j2) * length_cmap/length_fmap;
            end
        end
        

        % Dot product fields
        dmaps = zeros(11,9);
        dflowmap = zeros(11,9);
        dcellmap = zeros(11,9);
        dmaps(:,:)  = dot(flowmap(:,:,9:17), cmap(:,:,9:17));
        dflowmap(:,:)  = dot(flowmap(:,:,9:17), flowmap(:,:,9:17));
        dcellmap(:,:)  = dot(cmap(:,:,9:17), cmap(:,:,9:17));
        
        % Normalise dor products to sum for directions only
        if (directionscheme == 1)
            for i2 = 1:11
                for j2 = 1:9
                    if (dmaps(i2,j2) > 0)
                        dmaps(i2,j2) = dmaps(i2,j2)./sqrt(dflowmap(i2,j2).*dcellmap(i2,j2));
                    end
                end
            end
        end
        
        % Sum dot product fields along azimuths, into arrays 9 long
        azsumdmaps = sum(dmaps.* cos(el(:,9:17)*pi/180));

        % Sum azimuth sums
        sumdmaps = sum(azsumdmaps);

        % Cos corrected norm value
        dmaps_norm = sum(sum(ones(11,9).*cos(el(:,9:17)*pi/180)));
        
        % Normalise for number of points
        sumdmaps = sumdmaps / dmaps_norm;
        
        % Calculate score
        score(i,j) = sumdmaps;
        if (score(i,j) > maxscore)
            maxscore = score(i,j);
            maxaxis(1) = azi;
            maxaxis(2) = elj;
        end
    end
end

if (figureflag == 1)
    hold off
    quiver(az,el,cellmap{1},cellmap{2},'b')
    hold on
    if (rotation == 1)
        fmap = Generate_rotation_flow_map(maxaxis(1),maxaxis(2),0.1);
    else
        fmap = Generate_translation_flow_map(maxaxis(1),maxaxis(2),0.1);
    end
    fu(:,:) = fmap(1,:,:);
    fv(:,:) = fmap(2,:,:);
    quiver(az,el,fu,fv,'r')
end

% Attribute cell axis of rotation or translation
response{1} = maxaxis;
response{2} = maxscore;
response{3} = 0;
response{4} = score;
response{5} = cellmap{1};
response{6} = cellmap{2};
