function response = interpmap(mu,mv)
% assumes u,v is 1x66
% if not, checks in case it is already 11x17

x = zeros(66,1);
y = zeros(66,1);
x(1:7)   = [-120,-90,-45,0,45,90,120];
x(8:16)  = -120:30:120;
x(17:33) = -120:15:120;
x(34:50) = -120:15:120;
x(51:59) = -120:30:120;
x(60:66) = [-120,-90,-45,0,45,90,120];
y(1:7)   = 75;
y(8:16)  = 45;
y(17:33) = 15;
y(34:50) = -15;
y(51:59) = -45;
y(60:66) = -70;

[xi,yi] = meshgrid(-120:15:120,75:-30:-75);
yi(6,:) = -70;

if (length(mu) == 66)
    ui = zeros(6,17);
    vi = zeros(6,17);
    ui(1,:) = interp1(x(1:7),mu(1:7),(-120:15:120));
    vi(1,:) = interp1(x(1:7),mv(1:7),(-120:15:120));
    ui(2,:) = interp1(x(8:16),mu(8:16),(-120:15:120));
    vi(2,:) = interp1(x(8:16),mv(8:16),(-120:15:120));
    ui(3,:) = mu(17:33);
    vi(3,:) = mv(17:33);
    ui(4,:) = mu(34:50);
    vi(4,:) = mv(34:50);
    ui(5,:) = interp1(x(51:59),mu(51:59),(-120:15:120));
    vi(5,:) = interp1(x(51:59),mv(51:59),(-120:15:120));
    ui(6,:) = interp1(x(60:66),mu(60:66),(-120:15:120));
    vi(6,:) = interp1(x(60:66),mv(60:66),(-120:15:120));

    [xi2,yi2] = meshgrid(-120:15:120,75:-15:-75);
    yi2(11,:) = -70;
    ui2 = zeros(11,17);
    vi2 = zeros(11,17);
    ui2 = interp2(xi,yi,ui,xi2,yi2);
    vi2 = interp2(xi,yi,vi,xi2,yi2);
elseif (length(mu) == 17)
    ui2 = mu;
    vi2 = mv;
    % disp('already 11 x 17')
end

response{1} = ui2;
response{2} = vi2;
