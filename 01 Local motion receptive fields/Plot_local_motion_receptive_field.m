function Plot_local_motion_receptive_field(rec_id, cell_type)

if (cell_type == 1)
    % VT1 - wide field maps
    datdir{1} = '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-2/''';
    datdir{2} = '''../Data/Local motion receptive fields/VT1/Wide field maps/08 Jan 2014-1/''';
    datdir{3} = '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-1/''';
    datdir{4} = '''../Data/Local motion receptive fields/VT1/Wide field maps/10 Jan 2014-1/''';
    datdir{5} = '''../Data/Local motion receptive fields/VT1/Wide field maps/14 Jan 2014-1/''';
    datdir{6} = '''../Data/Local motion receptive fields/VT1/Wide field maps/21 Jan 2014-3/''';
    datdir{7} = '''../Data/Local motion receptive fields/VT1/Wide field maps/16 Apr 2008/''';
    datdir{8} = '''../Data/Local motion receptive fields/VT1/Wide field maps/29 Apr 2008/''';
    datdir{9} = '''../Data/Local motion receptive fields/VT1/Wide field maps/21 Mar 2007/''';
    datdir{10}= '''../Data/Local motion receptive fields/VT1/Wide field maps/25 Feb 2014-3/''';
    Ncells = length(datdir);
    s(1:Ncells) = 0.6;
    
elseif (cell_type == 2)
    % VT1 - local receptive field maps
    datdir{1} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/15 Jan 2014-1/''';
    datdir{2} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-2/''';
    datdir{3} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-4/''';
    datdir{4} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-7/''';
    datdir{5} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-1/''';
    datdir{6} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-4/''';
    datdir{7} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/11 Feb 2014-3/''';
    datdir{8} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/25 Feb 2014-4/''';
    datdir{9} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/27 Feb 2014-1/''';
    N_cells = length(datdir);
    s(1:N_cells) = 0.6;
    
elseif (cell_type == 3)
    % VT2
    datdir{1} = '''../Data/Local motion receptive fields/VT2/Wide field maps/04 Jun 2007-1/''';
    datdir{2} = '''../Data/Local motion receptive fields/VT2/Wide field maps/11 Dec 2007-1/''';
    datdir{3} = '''../Data/Local motion receptive fields/VT2/Wide field maps/26 Feb 2014-3/'''; % Incomplete
    Ncells = length(datdir);
    s(1:Ncells) = 0.6;
    
elseif (cell_type == 4)
    % VT3
    datdir{1} = '''../Data/Local motion receptive fields/VT3/Wide field maps/30 Apr 2008-1/''';
    datdir{2} = '''../Data/Local motion receptive fields/VT3/Wide field maps/12 Feb 2014-1/''';
    datdir{3} = '''../Data/Local motion receptive fields/VT3/Wide field maps/26 Feb 2014-1/''';
    datdir{4} = '''../Data/Local motion receptive fields/VT3/Wide field maps/11 Mar 2014-3/''';
    datdir{5} = '''../Data/Local motion receptive fields/VT3/Wide field maps/17 Mar 2014-4/''';
    Ncells = length(datdir);
    s(1:Ncells) = 0.6;
    
elseif (cell_type == 5)
    % Hu
    datdir{1} = '''../Data/Local motion receptive fields/Hu/Wide field maps/08 Apr 2008-1/''';
    datdir{2} = '''../Data/Local motion receptive fields/Hu/Wide field maps/09 Nov 2009-1/''';
    datdir{3} = '''../Data/Local motion receptive fields/Hu/Wide field maps/10 Nov 2009-1/''';
    datdir{4} = '''../Data/Local motion receptive fields/Hu/Wide field maps/19 Feb 2007-1/''';
    datdir{5} = '''../Data/Local motion receptive fields/Hu/Wide field maps/21 Feb 2007-1/''';
    datdir{6} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Feb 2007-2/''';
    datdir{7} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Mar 2007-1/''';
    datdir{8} = '''../Data/Local motion receptive fields/Hu/Wide field maps/28 Nov 2006-1/''';
    Ncells = length(datdir);
    s(1:Ncells) = 0.6;
    
end

% Load data
eval(['load ' datdir{rec_id} 'Map_data.mat']);

if (cell_type == 3 && rec_id == 1) % Flip round the map
    u2 = zeros(size(u));
    v2 = zeros(size(v));
    u2(1:7) = u(7:-1:1);
    u2(8:16) = u(16:-1:8);
    u2(17:33) = u(33:-1:17);
    u2(34:50) = u(50:-1:34);
    u2(51:59) = u(59:-1:51);
    u2(60:66) = u(66:-1:60);
    v2(1:7) = v(7:-1:1);
    v2(8:16) = v(16:-1:8);
    v2(17:33) = v(33:-1:17);
    v2(34:50) = v(50:-1:34);
    v2(51:59) = v(59:-1:51);
    v2(60:66) = v(66:-1:60);
    u = -u2;
    v = v2;
end


figure(rec_id + cell_type*100)
subplot(1,1,1)
quiver(x, y, u, v, s(rec_id));
set(gca,'tickdir','out')
set(gca,'xtick',-120:30:120);
set(gca,'ytick',[-70 -45 -15 15 45 75]);
xlabel(texlabel('azimuth ( ^o)'));
ylabel(texlabel('elevation ( ^o)'));
axis equal
axis([-140 140 -90 90])

