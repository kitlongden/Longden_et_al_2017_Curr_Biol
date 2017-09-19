function [az,el,daz,del] = Generate_translation_flow_map(azimuth,elevation,speed)

% azimuth = 0;
% elevation = 0;
% speed = 0.1;

% Specify axis of translation
axtr = [azimuth,elevation];
axtr(1) = -axtr(1) + 180;
axtr(2) = -axtr(2);

% Convert axis of translation to cartesian coordinates
axtr = axtr*2*pi/360;
xtr = cos(axtr(2)).*sin(axtr(1))*speed; 
ytr = cos(axtr(2)).*cos(axtr(1))*speed;
ztr = sin(axtr(2))*speed;

% Azimuth and elevation for interpolated maps
[az,el] = meshgrid(-120:15:120,75:-15:-75);
el(11,:) = -70;

% Conversion of azimuth and elevation to cartesian coordinates
% N.b. the azimuth runs counterclockwise in matlab
az2 = az*2*pi/360;
el2 = el*2*pi/360;
x = cos(el2).*sin(-az2); 
y = cos(el2).*cos(az2);
z = sin(el2);

% Translate all the azimuth elevation pairs in cartesian coordinates
xt =  x + xtr;
yt =  y + ytr;
zt =  z + ztr;

% Calculate the azimuth elevation angles of the translated coordinates
elt = atan2(zt, sqrt(xt.^2 + yt.^2));
elt = elt*360/(2*pi);
azt = -atan2(xt,yt);
azt = (azt*360/(2*pi));

daz = azt - az;
del = elt - el;

% quiver(az,el,daz,del,0,'k');

response = zeros(2,11,17);
response(1,:,:) = daz;
response(2,:,:) = del;
