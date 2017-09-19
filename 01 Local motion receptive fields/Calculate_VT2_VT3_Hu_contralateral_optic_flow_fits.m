function Calculate_VT2_VT3_Hu_contralateral_optic_flow_fits()

disp('VT3')
datdir{1} = '''../Data/Local motion receptive fields/VT3/Wide field maps/30 Apr 2008-1/''';
datdir{2} = '''../Data/Local motion receptive fields/VT3/Wide field maps/12 Feb 2014-1/''';
datdir{3} = '''../Data/Local motion receptive fields/VT3/Wide field maps/26 Feb 2014-1/''';
datdir{4} = '''../Data/Local motion receptive fields/VT3/Wide field maps/11 Mar 2014-3/''';
datdir{5} = '''../Data/Local motion receptive fields/VT3/Wide field maps/30 Apr 2008-1/''';
datadir = '../Data/''Local motion receptive fields''/VT3/Wide'' field maps''/';
fdir{1} = '30'' Apr 2008-1''/';
fdir{2} = '12'' Feb 2014-1''/';
fdir{3} = '26'' Feb 2014-1''/';
fdir{4} = '11'' Mar 2014-3''/';
fdir{5} = '30'' Apr 2008-1''/';
plot_fit_it(datadir,fdir,2);
clear datdir fdir datadir

disp('VT2')
datdir{1} = '''../Data/Local motion receptive fields/VT2/Wide field maps/04 Jun 2007-1/''';
datdir{2} = '''../Data/Local motion receptive fields/VT2/Wide field maps/11 Dec 2007-1/''';
datdir{3} = '''../Data/Local motion receptive fields/VT2/Wide field maps/26 Feb 2014-3/''';
datadir = '../Data/''Local motion receptive fields''/VT2/Wide'' field maps''/';
fdir{1} = '04'' Jun 2007-1''/';
fdir{2} = '11'' Dec 2007-1''/';
fdir{3} = '26'' Feb 2014-3''/';
plot_fit_it(datadir,fdir,1);
clear datdir fdir datadir

disp('Hu')
% Ventral visual field
datdir{1} = '''../Data/Local motion receptive fields/Hu/Wide field maps/08 Apr 2008-1/''';
datdir{2} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Feb 2007-2/''';
% Dorsal and ventral
datdir{3} = '''../Data/Local motion receptive fields/Hu/Wide field maps/19 Feb 2007-1/''';
datdir{4} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Nov 2006-1/''';
% Dorsal visual field
datdir{5} = '''../Data/Local motion receptive fields/Hu/Wide field maps/09 Nov 2009-1/''';
datdir{6} = '''../Data/Local motion receptive fields/Hu/Wide field maps/10 Nov 2009-1/''';
datdir{7} = '''../Data/Local motion receptive fields/Hu/Wide field maps/21 Feb 2007-1/''';
datdir{8} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Mar 2007-1/''';
datadir = '../Data/''Local motion receptive fields''/Hu/Wide'' field maps''/';
fdir{1} = '08'' Apr 2008-1''/';
fdir{2} = '28'' feb 2007-2''/';
fdir{3} = '19'' Feb 2007-1''/';
fdir{4} = '28'' Nov 2006-1''/';
fdir{5} = '09'' Nov 2009-1''/';
fdir{6} = '10'' Nov 2009-1''/';
fdir{7} = '21'' Feb 2007-1''/';
fdir{8} = '28'' Mar 2007-1''/';
plot_fit_it(datadir,fdir,3);
clear datdir fdir datadir

end


function plot_fit_it(datadir,fdir,cell_no)


normscheme = 3;
directionscheme = 1;
fname = 'Map_data.mat';

Ncells = length(fdir);

raxes = zeros(Ncells,2);
taxes = zeros(Ncells,2);
rscore = zeros(Ncells,1);
tscore = zeros(Ncells,1);
mu = zeros(Ncells,11,17);
mv = zeros(Ncells,11,17);
figureflag = 0;
for n = 1:Ncells
    disp(['Cell number ' num2str(n) '...'])
    eval(['load ' datadir fdir{n} fname]);
    
    r = Fit_rot_trans(1,normscheme,directionscheme,-180,-90,u,v,360,180,figureflag); % 110 125
    t = Fit_rot_trans(2,normscheme,directionscheme,-180,-90,u,v,360,180,figureflag); % 80 70
    
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

if (cell_no == 1)
    save ../Data/'Local motion receptive fields'/VT3/contralateral_optic_flow_fits_calculated.mat raxes taxes rscore tscore
elseif (cell_no == 2)
    save ../Data/'Local motion receptive fields'/VT2/contralateral_optic_flow_fits_calculated.mat raxes taxes rscore tscore
elseif (cell_no == 3)
    save ../Data/'Local motion receptive fields'/Hu/contralateral_optic_flow_fits_calculated.mat raxes taxes rscore tscore
end

end

