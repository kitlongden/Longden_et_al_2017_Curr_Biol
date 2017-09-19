function response = Generate_rotation_flow_map(azimuth,elevation,th)

% Specify axis of rotation and angle
axrot = [azimuth,elevation];


% Convert axis of rotation to cartesian coordinates
axrot = axrot*2*pi/360;
xrot = cos(axrot(2)).*cos(-axrot(1)); 
yrot = cos(axrot(2)).*sin(-axrot(1));
zrot = sin(axrot(2));

% Calculate rotation matrix
M = zeros(3,3);
M(1,1) = cos(th) + (1-cos(th))*(xrot^2);
M(1,2) = (1-cos(th))*xrot*yrot - sin(th)*zrot;
M(1,3) = (1-cos(th))*xrot*zrot + sin(th)*yrot;
M(2,1) = (1-cos(th))*yrot*xrot +sin(th)*zrot;
M(2,2) = cos(th) + (1-cos(th))*(yrot^2);
M(2,3) = (1-cos(th))*yrot*zrot - sin(th)*xrot;
M(3,1) = (1-cos(th))*zrot*xrot - sin(th)*yrot;
M(3,2) = (1-cos(th))*zrot*yrot + sin(th)*xrot;
M(3,3) = cos(th) + (1-cos(th))*(zrot^2);

% Azimuth and elevation for interpolated maps
[az,el] = meshgrid(-120:15:120,75:-15:-75);
el(11,:) = -70;

% Conversion of azimuth and elevation to cartesian coordinates
% N.b. the azimuth runs counterclockwise in matlab
az2 = az*2*pi/360;
el2 = el*2*pi/360;
x = cos(el2).*cos(-az2); 
y = cos(el2).*sin(-az2);
z = sin(el2);

% Rotate all the azimuth elevation pairs in cartesian coordinates
xr = zeros(11,17);
yr = zeros(11,17);
zr = zeros(11,17);
for i=1:11
    for j=1:17
        temp = M*[x(i,j);y(i,j);z(i,j)];
        xr(i,j) = temp(1);
        yr(i,j) = temp(2);
        zr(i,j) = temp(3);
    end
end

% Convert rotated coordinates back to azimuth elevation angles
elr = atan2(zr, sqrt(xr.^2 + yr.^2));
elr = elr*360/(2*pi);
azr = -atan2(yr,xr);
azr = azr*360/(2*pi);

response = zeros(2,11,17);
daz = azr-az;
del = elr-el;

% quiver(az,el,daz,del,0,'k');

response(1,:,:) = daz;
response(2,:,:) = del;

