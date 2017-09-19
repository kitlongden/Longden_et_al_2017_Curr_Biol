
% Parameters for <Fit_rot_trans_ipsi()>
normscheme = 3;
directionscheme = 1;

% Data
datadir = '../Data/''Local motion receptive fields''/VT1/''Wide field maps''/';
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
fdir{10}= '03'' Dec 2002''/';

Ncells = length(fdir);

raxes = zeros(Ncells,2);
taxes = zeros(Ncells,2);
rscore = zeros(Ncells,1);
tscore = zeros(Ncells,1);
mu = zeros(Ncells,11,17);
mv = zeros(Ncells,11,17);
figureflag = 0;

for n = 1:Ncells
    disp(n)
    if (n <= 9)
        eval(['load ' datadir fdir{n} fname]);
    else
        eval(['load ' datadir fdir{n} 'dot_mapping1_03-Dec-2002-SJH/Data_LPDS_xyuv.mat']);
    end
    r = Fit_rot_trans_ipsi(1,normscheme,directionscheme,-70,-55,u,v,110,125,figureflag); % 110 125
    t = Fit_rot_trans_ipsi(2,normscheme,directionscheme,  20,-10,u,v,80,70,figureflag); % 80 70
    raxes(n,:) = r{1};
    taxes(n,:) = t{1};
    rscore(n) = r{2};
    tscore(n) = t{2};
    mu(n,:,:) = r{5};
    mv(n,:,:) = r{6};
end

for n = 1:Ncells
    if (taxes(n) > 180)
        taxes(n) = taxes(n) - 360;
    end
end

save ../Data/'Local motion receptive fields'/VT1/ipsilateral_optic_flow_fits_calculated.mat raxes taxes rscore tscore