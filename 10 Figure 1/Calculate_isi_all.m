function Calculate_isi_all()

refract_period = 1.35;
samp_freq = 20;
burst_thr = 5;

calc_flag         = 1; % n=108 


%% Direction tuning
if (calc_flag == 1)
    disp('1 Direction tuning...')
    set_num = 1;
    N_cells = 6;
    N_trials = 20;
    N_directions = 16;
    datdir{1} = '''../Data/Direction tuning/15 Jan 2014-2/bindata/''';
    datdir{2} = '''../Data/Direction tuning/21 Jan 2014-1/bindata/''';
    datdir{3} = '''../Data/Direction tuning/21 Jan 2014-5/bindata/''';
    datdir{4} = '''../Data/Direction tuning/21 Jan 2014-6/bindata/''';
    datdir{5} = '''../Data/Direction tuning/22 Jan 2014-2/bindata/''';
    datdir{6} = '''../Data/Direction tuning/22 Jan 2014-3/bindata/''';
    calculate_isi(N_cells,N_trials,N_directions,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Size tuning
if (calc_flag == 1)
    disp('2 Size tuning...')
    set_num = 2;
    N_cells = 10;
    N_trials = 30;
    N_size = 20;
    datdir{1} = '''../Data/Size tuning/14 Apr 2014-2/bindata/''';
    datdir{2} = '''../Data/Size tuning/17 Mar 2014-2/bindata/''';
    datdir{3} = '''../Data/Size tuning/19 Mar 2014-2/bindata/''';
    datdir{4} = '''../Data/Size tuning/24 Mar 2014-2/bindata/''';
    datdir{5} = '''../Data/Size tuning/24 Mar 2014-3/bindata/''';
    datdir{6} = '''../Data/Size tuning/26 Mar 2014-3/bindata/''';
    datdir{7} = '''../Data/Size tuning/19 Mar 2014-6/bindata/''';
    datdir{8} = '''../Data/Size tuning/21 Mar 2014-2/bindata/''';
    datdir{9} = '''../Data/Size tuning/25 Mar 2014-1/bindata/''';
    datdir{10}= '''../Data/Size tuning/27 Mar 2014-2/bindata/''';
    calculate_isi(N_cells,N_trials,N_size,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end


%% Spatial frequency vs temporal frequency
if (calc_flag == 1)
    disp('3 Spatial frequency vs temporal frequency...')
    set_num = 3;
    N_cells = 10;
    N_trials = 25;
    N_tf = 8;
    N_sw = 6;
    datdir{1} =   '''../Data/Spatial wavelength and temporal frequency tuning/06 May 2014-3/bindata/''';
    datdir{2} =   '''../Data/Spatial wavelength and temporal frequency tuning/14 May 2014-7/bindata/''';
    datdir{3} =   '''../Data/Spatial wavelength and temporal frequency tuning/15 May 2014-2/bindata/''';
    datdir{4} =   '''../Data/Spatial wavelength and temporal frequency tuning/21 Apr 2014-3/bindata/''';
    datdir{5} =   '''../Data/Spatial wavelength and temporal frequency tuning/22 Apr 2014-3/bindata/''';
    datdir{6} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-10/bindata/''';
    datdir{7} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-2/bindata/''';
    datdir{8} =   '''../Data/Spatial wavelength and temporal frequency tuning/23 Apr 2014-7/bindata/''';
    datdir{9} =   '''../Data/Spatial wavelength and temporal frequency tuning/24 Apr 2014-2/bindata/''';
    datdir{10} =  '''../Data/Spatial wavelength and temporal frequency tuning/24 Apr 2014-5/bindata/''';
    calculate_isi(N_cells,N_trials,N_sw,N_tf,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Parallax motion
if (calc_flag == 1)
    disp('4 Parallax motion...')
    set_num = 4;
    N_cells = 10;
    N_trials = 10;
    N_tf = 9;
    datdir{1} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-1/bindata/''';
    datdir{2} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-2/bindata/''';
    datdir{3} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-4/bindata/''';
    datdir{4} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-6/bindata/''';
    datdir{5} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-8/bindata/''';
    datdir{6} = '''../Data/Motion in stimulus center and surround/VT1/28 Jan 2014-10/bindata/''';
    datdir{7} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-1/bindata/''';
    datdir{8} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-5/bindata/''';
    datdir{9} = '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-7/bindata/''';
    datdir{10}= '''../Data/Motion in stimulus center and surround/VT1/30 Jan 2014-9/bindata/''';
    calculate_isi(N_cells,N_trials,N_tf,N_tf,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Loom
if (calc_flag == 1)
    disp('5 Loom...')
    set_num = 5;
    N_cells = 11;
    N_trials = 30;
    N_looms = 8;
    datdir{1}  = '''../Data/Loom/17 Mar 2014-3/bindata/''';
    datdir{2}  = '''../Data/Loom/19 Mar 2014-3/bindata/''';
    datdir{3}  = '''../Data/Loom/26 Mar 2014-2/bindata/''';
    datdir{4}  = '''../Data/Loom/26 Mar 2014-4/bindata/''';
    datdir{5}  = '''../Data/Loom/27 Mar 2014-1/bindata/''';
    datdir{6}  = '''../Data/Loom/29 Mar 2014-1/bindata/''';
    datdir{7}  = '''../Data/Loom/29 Mar 2014-2/bindata/''';
    datdir{8}  = '''../Data/Loom/31 Mar 2014-1/bindata/''';
    datdir{9}  = '''../Data/Loom/31 Mar 2014-2/bindata/''';
    datdir{10} = '''../Data/Loom/31 Mar 2014-3/bindata/''';
    datdir{11} = '''../Data/Loom/14 Apr 2014-3/bindata/''';
    calculate_isi(N_cells,N_trials,N_looms,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Multiple depth planes 1
if (calc_flag == 1)
    disp('6 Multiple depth planes 1...')
    set_num = 6;
    N_cells = 9;
    N_trials = 30;
    N_step_heights = 8;
    datdir{1} = '''../Data/Multiple depth planes/Sine wave background/21 Apr 2014-2/bindata/''';
    datdir{2} = '''../Data/Multiple depth planes/Sine wave background/22 Apr 2014-2/bindata/''';
    datdir{3} = '''../Data/Multiple depth planes/Sine wave background/23 Apr 2014-4/bindata/''';
    datdir{4} = '''../Data/Multiple depth planes/Sine wave background/23 Apr 2014-6/bindata/''';
    datdir{5} = '''../Data/Multiple depth planes/Sine wave background/23 Apr 2014-9/bindata/''';
    datdir{6} = '''../Data/Multiple depth planes/Sine wave background/24 Apr 2014-1/bindata/''';
    datdir{7} = '''../Data/Multiple depth planes/Sine wave background/24 Apr 2014-4/bindata/''';
    datdir{8} = '''../Data/Multiple depth planes/Sine wave background/14 May 2014-5/bindata/''';
    datdir{9} = '''../Data/Multiple depth planes/Sine wave background/14 May 2014-6/bindata/''';
    calculate_isi(N_cells,N_trials,N_step_heights,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Multiple depth planes 2
if (calc_flag == 1)
    disp('7 Multiple depth planes 2...')
    set_num = 7;
    N_cells = 10;
    N_trials = 30;
    N_step_heights = 8;
    datdir{1}  = '''../Data/Multiple depth planes/Gray background/20 Apr 2014-7/bindata/''';
    datdir{2}  = '''../Data/Multiple depth planes/Gray background/21 Apr 2014-1/bindata/''';
    datdir{3}  = '''../Data/Multiple depth planes/Gray background/22 Apr 2014-1/bindata/''';
    datdir{4}  = '''../Data/Multiple depth planes/Gray background/23 Apr 2014-3/bindata/''';
    datdir{5}  = '''../Data/Multiple depth planes/Gray background/23 Apr 2014-5/bindata/''';
    datdir{6}  = '''../Data/Multiple depth planes/Gray background/23 Apr 2014-8/bindata/''';
    datdir{7}  = '''../Data/Multiple depth planes/Gray background/24 Apr 2014-3/bindata/''';
    datdir{8}  = '''../Data/Multiple depth planes/Gray background/06 May 2014-2/bindata/''';
    datdir{9}  = '''../Data/Multiple depth planes/Gray background/06 May 2014-5/bindata/''';
    datdir{10} = '''../Data/Multiple depth planes/Gray background/14 May 2014-2/bindata/''';
    calculate_isi(N_cells,N_trials,N_step_heights,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Flicker
if (calc_flag == 1)
    disp('8 Flicker...')
    set_num = 8;
    N_cells = 6;
    N_trials = 25;
    N_directions = 9;
    datdir{1} = '''../Data/Flicker/06 May 2014-1/bindata/''';
    datdir{2} = '''../Data/Flicker/06 May 2014-4/bindata/''';
    datdir{3} = '''../Data/Flicker/06 May 2014-6/bindata/''';
    datdir{4} = '''../Data/Flicker/14 May 2014-4/bindata/''';
    datdir{5} = '''../Data/Flicker/14 May 2014-8/bindata/''';
    datdir{6} = '''../Data/Flicker/15 May 2014-1/bindata/''';
    calculate_isi(N_cells,N_trials,N_directions,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% Maps
if (calc_flag == 1)
    disp('9 Maps...')
    set_num = 9;
    N_cells = 7;
    N_trials = 1;
    N_stim = 132;
    datdir{1} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-2/''';
    datdir{2} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/08 Jan 2014-1/''';
    datdir{3} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/09 Jan 2014-1/''';
    datdir{4} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/10 Jan 2014-1/''';
    datdir{5} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/14 Jan 2014-1/''';
    datdir{6} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/21 Jan 2014-3/''';
    datdir{7} =  '''../Data/Local motion receptive fields/VT1/Wide field maps/25 Feb 2014-3/''';
    calculate_isi(N_cells,N_trials,N_stim,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

%% White noiseflag == 1)
if (calc_flag == 1)
    disp('10 White noise...')
    set_num = 10;
    N_cells = 10;
    N_trials = 50;
    N_stim = 1;
    datdir{1} = '''../Data/Random motion stimulus/11 Feb 2014-5/bindata/''';
    datdir{2} = '''../Data/Random motion stimulus/27 Feb 2014-3/bindata/''';
    datdir{3} = '''../Data/Random motion stimulus/07 Mar 2014-1/bindata/''';
    datdir{4} = '''../Data/Random motion stimulus/07 Mar 2014-3/bindata/''';
    datdir{5} = '''../Data/Random motion stimulus/12 Mar 2014-1/bindata/''';
    datdir{6} = '''../Data/Random motion stimulus/12 Mar 2014-3/bindata/''';
    datdir{7} = '''../Data/Random motion stimulus/17 Mar 2014-1/bindata/''';
    datdir{8} = '''../Data/Random motion stimulus/19 Mar 2014-1/bindata/''';
    datdir{9} = '''../Data/Random motion stimulus/21 Mar 2014-1/bindata/''';
    datdir{10}= '''../Data/Random motion stimulus/24 Mar 2014-1/bindata/''';
    calculate_isi(N_cells,N_trials,N_stim,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end

if (calc_flag == 1)
    disp('11 RF Maps...')
    set_num = 11;
    N_cells = 9;
    N_trials = 1;
    N_stim = 22;
    datdir{1} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/15 Jan 2014-1/''';
    datdir{2} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-2/''';
    datdir{3} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-4/''';
    datdir{4} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/21 Jan 2014-7/''';
    datdir{5} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-1/''';
    datdir{6} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/22 Jan 2014-4/''';
    datdir{7} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/11 Feb 2014-3/''';
    datdir{8} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/25 Feb 2014-4/''';
    datdir{9} = '''../Data/Local motion receptive fields/VT1/Receptive field maps/27 Feb 2014-1/''';
    calculate_isi(N_cells,N_trials,N_stim,0,datdir,samp_freq,refract_period,burst_thr,set_num);
end
end

function calculate_isi(N_cells,N_trials,N_stim1,N_stim2,datdir,samp_freq,refract_period,burst_thr,set_num)

isi_all      = cell(N_cells,1);
isi_pre      = cell(N_cells,1);
isi_post     = cell(N_cells,1);
burst_num    = cell(N_cells,1);
burst2_preISI  = cell(N_cells,1);
burst2_postISI = cell(N_cells,1);
burst3_preISI  = cell(N_cells,1);
burst3_postISI = cell(N_cells,1);
burst_all_postISI = cell(N_cells,1);
burst_all_post2ISI = cell(N_cells,1);
for nc = 1:N_cells
    disp(['Cell...' num2str(nc)])
    
    isi_all{nc} = [];
    isi_pre{nc} = [];
    isi_post{nc} = [];
    burst_num{nc} = zeros(1,100);
    burst2_preISI{nc} = [];
    burst2_postISI{nc} = [];
    burst3_preISI{nc} = [];
    burst3_postISI{nc} = [];
    burst_all_postISI{nc} = [];
    burst_all_post2ISI{nc} = [];
    for nt = 1:N_trials
        % Load data
        if (set_num == 9)
            eval(['load ' datdir{nc} 'bindata.mat']);
            data = zeros(N_stim1,80000);
            data(1,:) = az120el75_cw_data_bin(:,1);
            data(2,:) = az120el75_ccw_data_bin(:,1);
            data(3,:) = az120el45_cw_data_bin(:,1);
            data(4,:) = az120el45_ccw_data_bin(:,1);
            data(5,:) = az120el15_cw_data_bin(:,1);
            data(6,:) = az120el15_ccw_data_bin(:,1);
            data(7,:) = az120eln15_cw_data_bin(:,1);
            data(8,:) = az120eln15_ccw_data_bin(:,1);
            data(9,:) = az120eln45_cw_data_bin(:,1);
            data(10,:) = az120eln45_ccw_data_bin(:,1);
            data(11,:) = az120eln70_cw_data_bin(:,1);
            data(12,:) = az120eln70_ccw_data_bin(:,1);
            data(13,:) = az105el15_cw_data_bin(:,1);
            data(14,:) = az105el15_ccw_data_bin(:,1);
            data(15,:) = az105eln15_cw_data_bin(:,1);
            data(16,:) = az105eln15_ccw_data_bin(:,1);
            data(17,:) = az90el75_cw_data_bin(:,1);
            data(18,:) = az90el75_ccw_data_bin(:,1);
            data(19,:) = az90el45_cw_data_bin(:,1);
            data(20,:) = az90el45_ccw_data_bin(:,1);
            data(21,:) = az90el15_cw_data_bin(:,1);
            data(22,:) = az90el15_ccw_data_bin(:,1);
            data(23,:) = az90eln15_cw_data_bin(:,1);
            data(24,:) = az90eln15_ccw_data_bin(:,1);
            data(25,:) = az90eln45_cw_data_bin(:,1);
            data(26,:) = az90eln45_ccw_data_bin(:,1);
            data(27,:) = az90eln70_cw_data_bin(:,1);
            data(28,:) = az90eln70_ccw_data_bin(:,1);
            data(29,:) = az75el15_cw_data_bin(:,1);
            data(30,:) = az75el15_ccw_data_bin(:,1);
            data(31,:) = az75eln15_cw_data_bin(:,1);
            data(32,:) = az75eln15_ccw_data_bin(:,1);
            data(33,:) = az60el45_cw_data_bin(:,1);
            data(34,:) = az60el45_ccw_data_bin(:,1);
            data(35,:) = az60el15_cw_data_bin(:,1);
            data(36,:) = az60el15_ccw_data_bin(:,1);
            data(37,:) = az60eln15_cw_data_bin(:,1);
            data(38,:) = az60eln15_ccw_data_bin(:,1);
            data(39,:) = az60eln45_cw_data_bin(:,1);
            data(40,:) = az60eln45_ccw_data_bin(:,1);
            data(41,:) = az45el75_cw_data_bin(:,1);
            data(42,:) = az45el75_ccw_data_bin(:,1);
            data(43,:) = az45el15_cw_data_bin(:,1);
            data(44,:) = az45el15_ccw_data_bin(:,1);
            data(45,:) = az45eln15_cw_data_bin(:,1);
            data(46,:) = az45eln15_ccw_data_bin(:,1);
            data(47,:) = az45eln70_cw_data_bin(:,1);
            data(48,:) = az45eln70_ccw_data_bin(:,1);
            data(49,:) = az30el15_cw_data_bin(:,1);
            data(50,:) = az30el15_ccw_data_bin(:,1);
            data(51,:) = az30el45_cw_data_bin(:,1);
            data(52,:) = az30el45_ccw_data_bin(:,1);
            data(53,:) = az30eln15_cw_data_bin(:,1);
            data(54,:) = az30eln15_ccw_data_bin(:,1);
            data(55,:) = az30eln45_cw_data_bin(:,1);
            data(56,:) = az30eln45_ccw_data_bin(:,1);
            data(57,:) = az15el15_cw_data_bin(:,1);
            data(58,:) = az15el15_ccw_data_bin(:,1);
            data(59,:) = az15eln15_cw_data_bin(:,1);
            data(60,:) = az15eln15_ccw_data_bin(:,1);
            data(61,:) = az0el75_cw_data_bin(:,1);
            data(62,:) = az0el75_ccw_data_bin(:,1);
            data(63,:) = az0el45_cw_data_bin(:,1);
            data(64,:) = az0el45_ccw_data_bin(:,1);
            data(65,:) = az0el15_cw_data_bin(:,1);
            data(66,:) = az0el15_ccw_data_bin(:,1);
            data(67,:) = az0eln15_cw_data_bin(:,1);
            data(68,:) = az0eln15_ccw_data_bin(:,1);
            data(69,:) = az0eln45_cw_data_bin(:,1);
            data(70,:) = az0eln45_ccw_data_bin(:,1);
            data(71,:) = az0eln70_cw_data_bin(:,1);
            data(72,:) = az0eln70_ccw_data_bin(:,1);
            data(73,:) = azn15el15_cw_data_bin(:,1);
            data(74,:) = azn15el15_ccw_data_bin(:,1);
            data(75,:) = azn15eln15_cw_data_bin(:,1);
            data(76,:) = azn15eln15_ccw_data_bin(:,1);
            data(77,:) = azn30el45_cw_data_bin(:,1);
            data(78,:) = azn30el45_ccw_data_bin(:,1);
            data(79,:) = azn30el15_cw_data_bin(:,1);
            data(80,:) = azn30el15_ccw_data_bin(:,1);
            data(81,:) = azn30eln15_cw_data_bin(:,1);
            data(82,:) = azn30eln15_ccw_data_bin(:,1);
            data(83,:) = azn30eln45_cw_data_bin(:,1);
            data(84,:) = azn30eln45_ccw_data_bin(:,1);
            data(85,:) = azn45el75_cw_data_bin(:,1);
            data(86,:) = azn45el75_ccw_data_bin(:,1);
            data(87,:) = azn45el15_cw_data_bin(:,1);
            data(88,:) = azn45el15_ccw_data_bin(:,1);
            data(89,:) = azn45eln15_cw_data_bin(:,1);
            data(90,:) = azn45eln15_ccw_data_bin(:,1);
            data(91,:) = azn45eln70_cw_data_bin(:,1);
            data(92,:) = azn45eln70_ccw_data_bin(:,1);
            data(93,:) = azn60el45_cw_data_bin(:,1);
            data(94,:) = azn60el45_ccw_data_bin(:,1);
            data(95,:) = azn60el15_cw_data_bin(:,1);
            data(96,:) = azn60el15_ccw_data_bin(:,1);
            data(97,:) = azn60eln15_cw_data_bin(:,1);
            data(98,:) = azn60eln15_ccw_data_bin(:,1);
            data(99,:) = azn60eln45_cw_data_bin(:,1);
            data(100,:) = azn60eln45_ccw_data_bin(:,1);
            data(101,:) = azn75el15_cw_data_bin(:,1);
            data(102,:) = azn75el15_ccw_data_bin(:,1);
            data(103,:) = azn75eln15_cw_data_bin(:,1);
            data(104,:) = azn75eln15_ccw_data_bin(:,1);
            data(105,:) = azn90el75_cw_data_bin(:,1);
            data(106,:) = azn90el75_ccw_data_bin(:,1);
            data(107,:) = azn90el45_cw_data_bin(:,1);
            data(108,:) = azn90el45_ccw_data_bin(:,1);
            data(109,:) = azn90el15_cw_data_bin(:,1);
            data(110,:) = azn90el15_ccw_data_bin(:,1);
            data(111,:) = azn90eln15_cw_data_bin(:,1);
            data(112,:) = azn90eln15_ccw_data_bin(:,1);
            data(113,:) = azn90eln45_cw_data_bin(:,1);
            data(114,:) = azn90eln45_ccw_data_bin(:,1);
            data(115,:) = azn90eln70_cw_data_bin(:,1);
            data(116,:) = azn90eln70_ccw_data_bin(:,1);
            data(117,:) = azn105el15_cw_data_bin(:,1);
            data(118,:) = azn105el15_ccw_data_bin(:,1);
            data(119,:) = azn105eln15_cw_data_bin(:,1);
            data(120,:) = azn105eln15_ccw_data_bin(:,1);
            data(121,:) = azn120el75_cw_data_bin(:,1);
            data(122,:) = azn120el75_ccw_data_bin(:,1);
            data(123,:) = azn120el45_cw_data_bin(:,1);
            data(124,:) = azn120el45_ccw_data_bin(:,1);
            data(125,:) = azn120el15_cw_data_bin(:,1);
            data(126,:) = azn120el15_ccw_data_bin(:,1);
            data(127,:) = azn120eln15_cw_data_bin(:,1);
            data(128,:) = azn120eln15_ccw_data_bin(:,1);
            data(129,:) = azn120eln45_cw_data_bin(:,1);
            data(130,:) = azn120eln45_ccw_data_bin(:,1);
            data(131,:) = azn120eln70_cw_data_bin(:,1);
            data(132,:) = azn120eln70_ccw_data_bin(:,1);
        elseif(set_num == 11)
            eval(['load ' datdir{nc} 'bindata.mat']);
            data = zeros(N_stim1,80000);
            data(1,:) = az15eln45_ccw_data_bin(:,1);
            data(2,:) = az45eln45_ccw_data_bin(:,1);
            data(3,:) = az75eln70_ccw_data_bin(:,1);
            data(4,:) = azn15eln70_ccw_data_bin(:,1);
            data(5,:) = az15eln45_cw_data_bin(:,1);
            data(6,:) = az45eln45_cw_data_bin(:,1);
            data(7,:) = az75eln70_cw_data_bin(:,1);
            data(8,:) = azn15eln70_cw_data_bin(:,1);
            data(9,:) = az105eln45_ccw_data_bin(:,1);
            data(10,:) = az15eln70_ccw_data_bin(:,1);
            data(11,:) = az60eln70_ccw_data_bin(:,1);
            data(12,:) = az105eln45_cw_data_bin(:,1);
            data(13,:) = az15eln70_cw_data_bin(:,1);
            data(14,:) = az60eln70_cw_data_bin(:,1);
            data(15,:) = az105eln70_ccw_data_bin(:,1);
            data(16,:) = az30eln70_ccw_data_bin(:,1);
            data(17,:) = az75eln45_ccw_data_bin(:,1);
            data(18,:) = azn15eln45_ccw_data_bin(:,1);
            data(19,:) = az105eln70_cw_data_bin(:,1);
            data(20,:) = az30eln70_cw_data_bin(:,1);
            data(21,:) = az75eln45_cw_data_bin(:,1);
            data(22,:) = azn15eln45_cw_data_bin(:,1);
        else
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat data;'])
        end
        if (N_stim2 == 0)
            for ns1 = 1:N_stim1
                % ISIs
                if (N_stim1 > 1)
                    isi = diff(find(data(ns1,:)>0))/samp_freq;
                else
                    isi = diff(find(data>0))'/samp_freq;
                end
                isi = isi(find(isi>refract_period));
                isi_all{nc} = cat(2,isi_all{nc},isi);
                if (length(isi)>2)
                    isi_pre{nc}  = cat(2,isi_pre{nc},isi(1:end-1));
                    isi_post{nc} = cat(2,isi_post{nc},isi(2:end));
                end
                % Burst number
                N_isi = length(isi);
                burst_count = 1;
                for ni = 1:N_isi
                    if (isi(ni) <=burst_thr)
                        burst_count = burst_count + 1;
                    else
                        burst_num{nc}(burst_count) = burst_num{nc}(burst_count) + 1;
                        % For 2 spike burst, burst pre and post ISI
                        if (burst_count == 2)
                            if (ni > 2)
                                burst2_preISI{nc}  = cat(2,burst2_preISI{nc},isi(ni-2));
                                burst2_postISI{nc} = cat(2,burst2_postISI{nc},isi(ni));
                            end
                        end
                        % For 3 spike burst, burst pre and post ISI
                        if (burst_count == 3)
                            if (ni > 3)
                                burst3_preISI{nc}  = cat(2,burst3_preISI{nc},isi(ni-3));
                                burst3_postISI{nc} = cat(2,burst3_postISI{nc},isi(ni));
                            end
                        end
                        % For all post burst ISI
                        burst_all_postISI{nc} = cat(2,burst_all_postISI{nc},isi(ni));
                        % For all post burst ISI
                        if (ni < N_isi)
                            burst_all_post2ISI{nc} = cat(2,burst_all_post2ISI{nc},isi(ni+1));
                        end
                        % Reset burst counter
                        burst_count = 1;
                    end
                end
            end
        else
            for ns1 = 1:N_stim1
                for ns2 = 1:N_stim2
                    % ISIs
                    isi = diff(find(data(ns1,ns2,:)>0))'/samp_freq;
                    isi = isi(find(isi>refract_period));
                    isi_all{nc} = cat(2,isi_all{nc},isi);
                    if (length(isi)>2)
                        isi_pre{nc}  = cat(2,isi_pre{nc},isi(1:end-1));
                        isi_post{nc} = cat(2,isi_post{nc},isi(2:end));
                    end
                    % Burst number
                    N_isi = length(isi);
                    burst_count = 1;
                    for ni = 1:N_isi
                        if (isi(ni) <=burst_thr)
                            burst_count = burst_count + 1;
                        else
                            burst_num{nc}(burst_count) = burst_num{nc}(burst_count) + 1;
                            % For 2 spike burst, burst pre and post ISI
                            if (burst_count == 2)
                                if (ni > 2)
                                    burst2_preISI{nc}  = cat(2,burst2_preISI{nc},isi(ni-2));
                                    burst2_postISI{nc} = cat(2,burst2_postISI{nc},isi(ni));
                                end
                            end
                            % For 3 spike burst, burst pre and post ISI
                            if (burst_count == 3)
                                if (ni > 3)
                                    burst3_preISI{nc}  = cat(2,burst3_preISI{nc},isi(ni-3));
                                    burst3_postISI{nc} = cat(2,burst3_postISI{nc},isi(ni));
                                end
                            end
                            % For all post burst ISI
                            burst_all_postISI{nc} = cat(2,burst_all_postISI{nc},isi(ni));
                            % For all post burst ISI
                            if (ni < N_isi)
                                burst_all_post2ISI{nc} = cat(2,burst_all_post2ISI{nc},isi(ni+1));
                            end
                            % Reset burst counter
                            burst_count = 1;
                        end
                    end
                end
            end
        end
    end
end

% Save
eval(['save ISI_data_' num2str(set_num) ' isi_all isi_pre isi_post burst_num burst2_preISI burst2_postISI burst3_preISI burst3_postISI burst_all_postISI burst_all_post2ISI'])
end

