% Stimulis are in units of pixel shifts. 25 Hz corresponds to a shift of
% 16.61 pixels in one frame.
tf_fact = 25/16.61;

first_flag = 1;
stim_flag  = 1;
mean_flag  = 1;
graph_flag = 1;

N_cells = 10;
N_trials = 50;
N_samples = 600000;
dt = 3000; % 150 ms

% The data format is as follows:
% data(sizes, samples).

datadir{1} = '''../Data/Random motion stimulus/11 Feb 2014-5/bindata/''';
datadir{2} = '''../Data/Random motion stimulus/27 Feb 2014-3/bindata/''';
datadir{3} = '''../Data/Random motion stimulus/07 Mar 2014-1/bindata/''';
datadir{4} = '''../Data/Random motion stimulus/07 Mar 2014-3/bindata/''';
datadir{5} = '''../Data/Random motion stimulus/12 Mar 2014-1/bindata/''';
datadir{6} = '''../Data/Random motion stimulus/12 Mar 2014-3/bindata/''';
datadir{7} = '''../Data/Random motion stimulus/17 Mar 2014-1/bindata/''';
datadir{8} = '''../Data/Random motion stimulus/19 Mar 2014-1/bindata/''';
datadir{9} = '''../Data/Random motion stimulus/21 Mar 2014-1/bindata/''';
datadir{10}= '''../Data/Random motion stimulus/24 Mar 2014-1/bindata/''';

if (first_flag == 1)
    disp('   Burst spike array...')
    burst_array = cell(N_cells,1);
    burst_array2 = cell(N_cells,1);
    fs = cell(N_cells,1);
    
    for nc = 1:N_cells
        disp(['Cell...' num2str(nc)])
        
        disp('   Preparing stimulus...')
        eval(['load ' datadir{nc} 'frameshift.mat;'])
        fs{nc} = zeros(N_trials,N_samples);
        for nt = 1:N_trials
            for i = 1:100
                fs{nc}(nt,i:100:N_samples) = frameshift(nt,:);
            end
        end
        
        disp('   Loading spike times...')
        id = cell(N_trials,1);
        for nt = 1:N_trials
            eval(['load ' datadir{nc} 'trial_' num2str(nt) '.mat;'])
            % The data format is: data(number of samples,1).
            id{nt} = (find(data(:) > 0))';
            id{nt} = id{nt} / 20; % conversion to ms
        end
        
        disp('   Identifying bursts...') % Labels spikes with # spikes in burst
        burst_array{nc} = classify_spikes_fn_STA(id, N_trials, N_samples);
        
        disp('   Eliminate all but the first spikes of bursts...')
        for nt = 1:N_trials
            data = burst_array{nc}(nt,:);
            id = find(data>0);
            N_id = length(id);
            for ni = 1:N_id
                burst_size = data(id(ni));
                if  (burst_size > 1)
                    data(id(ni+1:ni+burst_size-1)) = 0;
                end
            end
            burst_array2{nc}(nt,:) = data;
        end
    end
    
    % disp('    save r_burst data...')
    save -v7.3 sta_data burst_array burst_array2 fs
end    

if (stim_flag == 1)
    disp('   Stimulus arrays...')
    
    load sta_data burst_array2 fs
    stim1   = cell(N_cells,1);
    
    for nc = 1:N_cells
        disp(nc)
        disp('...Single spikes...')
        stim1{nc} = [];
        for nt = 1:N_trials
            data = burst_array2{nc}(nt,:);
            id  = find(data == 1);
            Nid = length(id);
            stim = zeros(Nid,dt);
            for i = 1:Nid
                if (id(i) > dt)
                    stim(i,:) = fs{nc}(nt,id(i)-dt+1:id(i));
                end
            end
            stim1{nc} = cat(1,stim1{nc},stim);
        end
        disp('...2 spike bursts...')
        stim2{nc} = [];
        for nt = 1:N_trials
            data = burst_array2{nc}(nt,:);
            id  = find(data == 2);
            Nid = length(id);
            stim = zeros(Nid,dt);
            for i = 1:Nid
                if (id(i) > dt)
                    stim(i,:) = fs{nc}(nt,id(i)-dt+1:id(i));
                end
            end
            stim2{nc} = cat(1,stim2{nc},stim);
        end
        disp('...3 spikes...')
        stim3{nc} = [];
        for nt = 1:N_trials
            data = burst_array2{nc}(nt,:);
            id  = find(data == 3);
            Nid = length(id);
            stim = zeros(Nid,dt);
            for i = 1:Nid
                if (id(i) > dt)
                    stim(i,:) = fs{nc}(nt,id(i)-dt+1:id(i));
                end
            end
            stim3{nc} = cat(1,stim3{nc},stim);
        end
        disp('...3+ spikes...')
        stim3p{nc} = [];
        for nt = 1:N_trials
            data = burst_array2{nc}(nt,:);
            id  = find(data >= 3);
            Nid = length(id);
            stim = zeros(Nid,dt);
            for i = 1:Nid
                if (id(i) > dt)
                    stim(i,:) = fs{nc}(nt,id(i)-dt+1:id(i));
                end
            end
            stim3p{nc} = cat(1,stim3p{nc},stim);
        end
        disp('...4 spikes...')
        stim4{nc} = [];
        for nt = 1:N_trials
            data = burst_array2{nc}(nt,:);
            id  = find(data == 4);
            Nid = length(id);
            stim = zeros(Nid,dt);
            for i = 1:Nid
                if (id(i) > dt)
                    stim(i,:) = fs{nc}(nt,id(i)-dt+1:id(i));
                end
            end
            stim4{nc} = cat(1,stim4{nc},stim);
        end
    end
   
    % Save 
    save -append sta_data stim1 stim2 stim3 stim3p stim4
end

    
if (mean_flag == 1)
    disp('   Mean stimulus arrays...')
    load sta_data stim1 stim2 stim3 stim3p stim4
    
    mstim1  = zeros(N_cells,dt);
    mstim2  = zeros(N_cells,dt);
    mstim3  = zeros(N_cells,dt);
    mstim3p = zeros(N_cells,dt);
    mstim4  = zeros(N_cells,dt);
    
    for nc = 1:N_cells
        mstim1(nc,:)  = mean(stim1{nc})*tf_fact;
        mstim2(nc,:)  = mean(stim2{nc})*tf_fact;
        mstim3(nc,:)  = mean(stim3{nc})*tf_fact;
        mstim3p(nc,:) = mean(stim3p{nc})*tf_fact;
        mstim4(nc,:)  = mean(stim4{nc})*tf_fact;
    end
    
    % Save 
    save -append sta_data mstim1 mstim2 mstim3 mstim3p mstim4
end
    
    
if (graph_flag == 1)
    disp('   Plot mean STAs.')
    
    load sta_data mstim1 mstim2 mstim3 mstim3p mstim4
    mn_mstim_all = zeros(5,dt);
    
    mn_mstim_all(1,:) = mean(mstim1);
    mn_mstim_all(2,:) = mean(mstim2);
    mn_mstim_all(3,:) = mean(mstim3);
    mn_mstim_all(4,:) = mean(mstim3p);
    mn_mstim_all(5,:) = mean(mstim4);
    
    se_mstim_all(1,:) = std(mstim1)/sqrt(N_cells);
    se_mstim_all(2,:) = std(mstim2)/sqrt(N_cells);
    se_mstim_all(3,:) = std(mstim3)/sqrt(N_cells);
    se_mstim_all(4,:) = std(mstim3p)/sqrt(N_cells);
    se_mstim_all(5,:) = std(mstim4)/sqrt(N_cells);
    
    figure(3)
    subplot(2,2,1)
    plot(-mn_mstim_all([1 2 3 5],:)')
    
    subplot(2,2,2)
    plot(-mn_mstim_all([1 2 4],:)')
    
    subplot(2,2,3)
    hold off
    errorbar((1:dt)/20,-mn_mstim_all(1,:)',-se_mstim_all(1,:)')
    hold on
    errorbar((1:dt)/20,-mn_mstim_all(2,:)',-se_mstim_all(2,:)')
    errorbar((1:dt)/20,-mn_mstim_all(3,:)',-se_mstim_all(3,:)')
    errorbar((1:dt)/20,-mn_mstim_all(5,:)',-se_mstim_all(5,:)')
    
    subplot(2,2,4)
    hold off
    errorbar((-dt+1:0)/20,-mn_mstim_all(4,:)',-se_mstim_all(4,:)')
    hold on
    errorbar((-dt+1:0)/20,-mn_mstim_all(2,:)',-se_mstim_all(2,:)')
    errorbar((-dt+1:0)/20,-mn_mstim_all(1,:)',-se_mstim_all(1,:)')
    
    
    
    % Save
    save -append sta_data mn_mstim_all se_mstim_all
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    








