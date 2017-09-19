
%% 29th October 2016

calc_flag = 1;

N_cells = 11;
N_trials = 30;
N_looms = 8;
N_samples = 80000;
N_ms = 8000;
samp_freq = 20000;
    
% in mm/s 
loom_speeds = [30 150 500 5000];

% 
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

if (calc_flag == 1)
    
    samp_1 = 1;
    samp_2 = 80000;
    samp_duration = (samp_2-samp_1)/samp_freq;
    stage1_1 = 1;
    stage1_2 = 40000;
    stage1_duration = (stage1_2-stage1_1)/samp_freq;
    stage2_1 = 40001;
    stage2_2 = 60000;
    stage2_duration = (stage2_2-stage2_1)/samp_freq;
    stage3_1 = 6001;
    stage3_2 = 70000;
    stage3_duration = (stage3_2-stage3_1)/samp_freq;
    stage4_1 = 70001;
    stage4_2 = 80000;
    stage4_duration = (stage4_2-stage4_1)/samp_freq;
    
    r_burst_1  = zeros(N_cells, N_looms, N_trials);
    r_burst_2  = zeros(N_cells, N_looms, N_trials);
    r_burst_3  = zeros(N_cells, N_looms, N_trials);
    r_stage1_1 = zeros(N_cells, N_looms, N_trials);
    r_stage1_2 = zeros(N_cells, N_looms, N_trials);
    r_stage1_3 = zeros(N_cells, N_looms, N_trials);
    r_stage2_1 = zeros(N_cells, N_looms, N_trials);
    r_stage2_2 = zeros(N_cells, N_looms, N_trials);
    r_stage2_3 = zeros(N_cells, N_looms, N_trials);
    r_stage3_1 = zeros(N_cells, N_looms, N_trials);
    r_stage3_2 = zeros(N_cells, N_looms, N_trials);
    r_stage3_3 = zeros(N_cells, N_looms, N_trials);
    r_stage4_1 = zeros(N_cells, N_looms, N_trials);
    r_stage4_2 = zeros(N_cells, N_looms, N_trials);
    r_stage4_3 = zeros(N_cells, N_looms, N_trials);
    
    burst_counts = zeros(N_cells, 10);
    barray_1 = cell(N_cells,1);
    barray_2 = cell(N_cells,2);
    barray_3 = cell(N_cells,3);
    
    gw = gausswin(800,2.5);
    gw = gw / sum(gw);
    fbarray_1 = cell(N_cells,1);
    fbarray_2 = cell(N_cells,2);
    fbarray_3 = cell(N_cells,3);
    
    afbarray_1 = zeros(N_looms,N_samples);
    afbarray_2 = zeros(N_looms,N_samples);
    afbarray_3 = zeros(N_looms,N_samples);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Loading spike times...')
        id = cell(N_looms, N_trials);
        for nt = 1:N_trials
            eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is as follows:
            % data(N_looms, N_samples).
            for nl = 1:N_looms
                id{nl,nt} = (find(data(nl,:) > 0))';
                id{nl,nt} = id{nl,nt} / 20; % conversion to ms
            end
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array = classify_spikes_fn_loom(id, N_looms, N_trials, N_samples);
        
        disp('   Calculate rate responses for burst categories...')
        for nl = 1:N_looms
            for nt = 1:N_trials
                % Total response
                r_burst_1(nc,nl,nt) = length(find(burst_array(nl,nt,samp_1:samp_2) == 1))/samp_duration;
                r_burst_2(nc,nl,nt) = length(find(burst_array(nl,nt,samp_1:samp_2) == 2))/samp_duration;
                r_burst_3(nc,nl,nt) = length(find(burst_array(nl,nt,samp_1:samp_2) >= 3))/samp_duration;
                % Stage1 response
                r_stage1_1(nc,nl,nt) = length(find(burst_array(nl,nt,stage1_1:stage1_2) == 1))/stage1_duration;
                r_stage1_2(nc,nl,nt) = length(find(burst_array(nl,nt,stage1_1:stage1_2) == 2))/stage1_duration;
                r_stage1_3(nc,nl,nt) = length(find(burst_array(nl,nt,stage1_1:stage1_2) >= 3))/stage1_duration;
                 % Stage2 response
                r_stage2_1(nc,nl,nt) = length(find(burst_array(nl,nt,stage2_1:stage2_2) == 1))/stage2_duration;
                r_stage2_2(nc,nl,nt) = length(find(burst_array(nl,nt,stage2_1:stage2_2) == 2))/stage2_duration;
                r_stage2_3(nc,nl,nt) = length(find(burst_array(nl,nt,stage2_1:stage2_2) >= 3))/stage2_duration;
                 % Stage3 response
                r_stage3_1(nc,nl,nt) = length(find(burst_array(nl,nt,stage3_1:stage3_2) == 1))/stage3_duration;
                r_stage3_2(nc,nl,nt) = length(find(burst_array(nl,nt,stage3_1:stage3_2) == 2))/stage3_duration;
                r_stage3_3(nc,nl,nt) = length(find(burst_array(nl,nt,stage3_1:stage3_2) >= 3))/stage3_duration;
                 % Stage4 response
                r_stage4_1(nc,nl,nt) = length(find(burst_array(nl,nt,stage4_1:stage4_2) == 1))/stage4_duration;
                r_stage4_2(nc,nl,nt) = length(find(burst_array(nl,nt,stage4_1:stage4_2) == 2))/stage4_duration;
                r_stage4_3(nc,nl,nt) = length(find(burst_array(nl,nt,stage4_1:stage4_2) >= 3))/stage4_duration;
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
        for nl = 1:N_looms
            fbarray_1{nc}(nl,:) = filtfilt(gw,1,barray_1{nc}(nl,:));
            fbarray_2{nc}(nl,:) = filtfilt(gw,1,barray_2{nc}(nl,:));
            fbarray_3{nc}(nl,:) = filtfilt(gw,1,barray_3{nc}(nl,:));
            afbarray_1(nl,:) = afbarray_1(nl,:) + fbarray_1{nc}(nl,:);
            afbarray_2(nl,:) = afbarray_2(nl,:) + fbarray_2{nc}(nl,:);
            afbarray_3(nl,:) = afbarray_3(nl,:) + fbarray_3{nc}(nl,:); 
        end
    end
    afbarray_1 = afbarray_1 / N_cells;
    afbarray_2 = afbarray_2 / N_cells;
    afbarray_3 = afbarray_3 / N_cells;
    
    disp('    save r_burst data...')
    save r_burst_data r_burst* burst_counts barray_* fbarray_* afbarray_* r_stage*
end