

% Variables
N_data = 132;
sample_rate = 20000; % Hz
duration = 4; % seconds
ms_factor = 1000;% miliseconds
burst_isi_criterion = 5; % ms; if ISI is less, it's a burst
N_max_spikes_per_burst = 10;

% Bin ISI variables
dbin = .5; % ms
N_bins = 2000;
for nb = 1:N_bins
    bin_min(nb) = (nb-1) * dbin;
    bin_max(nb) = nb * dbin;
end
bin_centres = ((0:dbin:(N_bins-1)*dbin)) + dbin/2;

% Cells
clear datdir
datdir{1} =  '''Data/VT1/Wide field maps/09 Jan 2014-2/''';
datdir{2} =  '''Data/VT1/Wide field maps/08 Jan 2014-1/''';
datdir{3} =  '''Data/VT1/Wide field maps/09 Jan 2014-1/''';
datdir{4} =  '''Data/VT1/Wide field maps/10 Jan 2014-1/''';
datdir{5} =  '''Data/VT1/Wide field maps/14 Jan 2014-1/''';
datdir{6} =  '''Data/VT1/Wide field maps/21 Jan 2014-3/''';
datdir{7} =  '''Data/VT1/Wide field maps/25 Feb 2014-3/''';

N_cells = length(datdir);
% N_cells = 1;
spike_burst_num = cell(N_cells,1);
bin_spike_burst_num = NaN(N_cells, N_max_spikes_per_burst);

% Array for Binning ISIs
bin_count = zeros(N_cells, N_bins);
norm_bin_count = zeros(N_cells, N_bins);

% Load data
for nc = 1:N_cells
    eval(['load ' datdir{nc} 'bindata.mat']);
    disp(['Cell...' num2str(nc)])
    
    data = zeros(N_data,sample_rate*duration);
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
    
    % Clear data at very start and end of trials that could have been
    % mid-burst.
    data(:,1:1+(burst_isi_criterion*20)) = 0;
    data(:,end-(burst_isi_criterion*20):end) = 0;
    
    % Extract ISIs
    isi = [];
    spike_burst_num{nc} = [];
    for nd = 1:N_data
        if (sum(data(nd,:)) > 0)
            spike_t = find(data(nd,:) > 0)/(sample_rate/ms_factor);
            N_spike_t = length(spike_t);
            % Pick out ISIs
            if (N_spike_t > 1)
                isi = [isi diff(spike_t)]; %#ok
            end
            % Pick out the number of spikes per burst
            if (N_spike_t == 1)
                spike_burst_num{nc} = [spike_burst_num{nc} 1]; %#ok
            else % Identify multiple spike bursts
                burst_count = 1;
                for ns = 1:N_spike_t-1
                    if (spike_t(ns+1) < spike_t(ns) + burst_isi_criterion)
                        burst_count = burst_count + 1;
                    else
                        spike_burst_num{nc} = [spike_burst_num{nc} burst_count]; %#ok
                        burst_count = 1;
                    end
                end
            end
        end
    end
    
    % Bin ISIs
    for nb = 1:N_bins
        bin_count(nc,nb) = length(find(isi >= bin_min(nb) & isi < bin_max(nb)));
    end
    norm_bin_count(nc,:) = bin_count(nc,:) / sum(bin_count(nc,:));
    
    % Bin spike burst numbers
    for nsb = 1:N_max_spikes_per_burst
        bin_spike_burst_num(nc,nsb) = length(find(spike_burst_num{nc} == nsb));
    end
    bin_spike_burst_num(nc,:) = bin_spike_burst_num(nc,:)/sum(bin_spike_burst_num(nc,:));
end

% Save isi data
save Data/VT1/map_isi_data.mat bin_count norm_bin_count bin_spike_burst_num spike_burst_num N_cells isi
