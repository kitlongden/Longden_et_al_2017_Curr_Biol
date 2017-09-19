
calc_flag = 1;

N_trials = 30;
N_step_heights = 8;
N_samples = 80000;
N_ms = 4000;
samp_freq = 20000;

% The data format is as follows:
% data(step_heights, samples).
clear datdir

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
N_cells = length(datdir);

if (calc_flag == 1)
    
    samp_1 = 1;
    samp_2 = 80000;
    samp_duration = (samp_2-samp_1)/samp_freq;
    stage1_1 = 1;
    stage1_2 = 40000;
    stage1_duration = (stage1_2-stage1_1)/samp_freq;
    stage2_1 = 40001;
    stage2_2 = 50000;
    stage2_duration = (stage2_2-stage2_1)/samp_freq;
    stage3_1 = 4001;
    stage3_2 = 60000;
    stage3_duration = (stage3_2-stage3_1)/samp_freq;
    stage4_1 = 60001;
    stage4_2 = 70000;
    stage4_duration = (stage4_2-stage4_1)/samp_freq;
    
    r_burst_1 = zeros(N_cells, N_step_heights, N_trials);
    r_burst_2 = zeros(N_cells, N_step_heights, N_trials);
    r_burst_3 = zeros(N_cells, N_step_heights, N_trials);
    r_stage1_1 = zeros(N_cells, N_step_heights, N_trials);
    r_stage1_2 = zeros(N_cells, N_step_heights, N_trials);
    r_stage1_3 = zeros(N_cells, N_step_heights, N_trials);
    r_stage2_1 = zeros(N_cells, N_step_heights, N_trials);
    r_stage2_2 = zeros(N_cells, N_step_heights, N_trials);
    r_stage2_3 = zeros(N_cells, N_step_heights, N_trials);
    r_stage3_1 = zeros(N_cells, N_step_heights, N_trials);
    r_stage3_2 = zeros(N_cells, N_step_heights, N_trials);
    r_stage3_3 = zeros(N_cells, N_step_heights, N_trials);
    r_stage4_1 = zeros(N_cells, N_step_heights, N_trials);
    r_stage4_2 = zeros(N_cells, N_step_heights, N_trials);
    r_stage4_3 = zeros(N_cells, N_step_heights, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    gw = gausswin(800,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_step_heights,N_samples);
    afbarray_2 = zeros(N_step_heights,N_samples);
    afbarray_3 = zeros(N_step_heights,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_step_heights, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(N_step_heights, N_samples).
            for nsh = 1:N_step_heights
                id{nsh,nt} = (find(data(nsh,:) > 0))';
                id{nsh,nt} = id{nsh,nt} / 20; % conversion to ms
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_depth(id, N_step_heights, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for nsh = 1:N_step_heights
            for nt = 1:N_trials
                % Total response
                r_burst_1(nc,nsh,nt) = length(find(burst_array(nsh,nt,samp_1:samp_2) == 1))/samp_duration;
                r_burst_2(nc,nsh,nt) = length(find(burst_array(nsh,nt,samp_1:samp_2) == 2))/samp_duration;
                r_burst_3(nc,nsh,nt) = length(find(burst_array(nsh,nt,samp_1:samp_2) >= 3))/samp_duration;
                % Stage1 response
                r_stage1_1(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage1_1:stage1_2) == 1))/stage1_duration;
                r_stage1_2(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage1_1:stage1_2) == 2))/stage1_duration;
                r_stage1_3(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage1_1:stage1_2) >= 3))/stage1_duration;
                 % Stage2 response
                r_stage2_1(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage2_1:stage2_2) == 1))/stage2_duration;
                r_stage2_2(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage2_1:stage2_2) == 2))/stage2_duration;
                r_stage2_3(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage2_1:stage2_2) >= 3))/stage2_duration;
                 % Stage3 response
                r_stage3_1(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage3_1:stage3_2) == 1))/stage3_duration;
                r_stage3_2(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage3_1:stage3_2) == 2))/stage3_duration;
                r_stage3_3(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage3_1:stage3_2) >= 3))/stage3_duration;
                 % Stage4 response
                r_stage4_1(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage4_1:stage4_2) == 1))/stage4_duration;
                r_stage4_2(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage4_1:stage4_2) == 2))/stage4_duration;
                r_stage4_3(nc,nsh,nt) = length(find(burst_array(nsh,nt,stage4_1:stage4_2) >= 3))/stage4_duration;
            end
        end
        
        disp('   Calculate ISI distributions...')
        for nsb = 1:10
            burst_counts(nc,nsb) = length(find(burst_array == nsb));
        end
        
        disp('   Calculate mean response traces...')
        burst_array_1 = (burst_array == 1)*samp_freq;
        burst_array_2 = (burst_array == 2)*samp_freq;
        burst_array_3 = (burst_array >= 3)*samp_freq;
        barray_1{nc} = squeeze(mean(burst_array_1,2));
        barray_2{nc} = squeeze(mean(burst_array_2,2));
        barray_3{nc} = squeeze(mean(burst_array_3,2));
        
        disp('   Filter mean response traces...')
        for nsh = 1:N_step_heights
            fbarray_1{nc}(nsh,:) = filtfilt(gw,1,barray_1{nc}(nsh,:));
            fbarray_2{nc}(nsh,:) = filtfilt(gw,1,barray_2{nc}(nsh,:));
            fbarray_3{nc}(nsh,:) = filtfilt(gw,1,barray_3{nc}(nsh,:));
            afbarray_1(nsh,:) = afbarray_1(nsh,:) + fbarray_1{nc}(nsh,:);
            afbarray_2(nsh,:) = afbarray_2(nsh,:) + fbarray_2{nc}(nsh,:);
            afbarray_3(nsh,:) = afbarray_3(nsh,:) + fbarray_3{nc}(nsh,:); 
        end
    end
    afbarray_1 = afbarray_1 / N_cells;
    afbarray_2 = afbarray_2 / N_cells;
    afbarray_3 = afbarray_3 / N_cells;
    
    disp('    save r_burst data...')
    save r_burst_data_grey r_burst* burst_counts barray_* fbarray_* afbarray_* r_stage*
end