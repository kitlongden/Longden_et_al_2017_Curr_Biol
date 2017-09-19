function Calculate_isi_all_V1()

%%
refract_period = 1.35;
samp_freq = 20;
burst_thr = 5;

N_trials = 10;
N_tf = 9;
datdir{1} = '''../Data/Motion in stimulus center and surround/V1/03 Feb 2014-1/bindata/''';
datdir{2} = '''../Data/Motion in stimulus center and surround/V1/03 Feb 2014-6/bindata/''';
datdir{3} = '''../Data/Motion in stimulus center and surround/V1/03 Feb 2014-8/bindata/''';
datdir{4} = '''../Data/Motion in stimulus center and surround/V1/04 Feb 2014-1/bindata/''';
datdir{5} = '''../Data/Motion in stimulus center and surround/V1/30 Jan 2014-3/bindata/''';
datdir{6} = '''../Data/Motion in stimulus center and surround/V1/31 Jan 2014-2/bindata/''';
N_cells = 6;
calculate_isi(N_cells,N_trials,N_tf,N_tf,datdir,samp_freq,refract_period,burst_thr);


%% Numbers of spikes per burst
burst_num_all_norm = [];
load data_V1.mat burst_num
N_cells = length(burst_num);
for nc = 1:N_cells
    if (sum(burst_num{nc}) > 0)
        burst_num_all_norm = cat(1,burst_num_all_norm,burst_num{nc}/sum(burst_num{nc}));
    end
end
save data_V1_all burst_num_all_norm;


%% Bin ISIs
disp('Binning ISIs...')
bin_max = 1000;
db = 0.5;
bins = 0:db:bin_max;
bins2 = db/2:db:bin_max-db/2;
N_bins = length(bins) - 1;
isi_norm = [];

load data_V1.mat
N_cells = length(isi_all);
B_bins = zeros(N_cells,N_bins);
for nc = 1:N_cells
    for nb = 1:N_bins
        id = find(isi_all{nc} > bins(nb) & isi_all{nc} <= bins(nb+1));
        B_bins(nc,nb) = length(id);
    end
    B_bins(nc,:) = B_bins(nc,:)/sum(B_bins(nc,:));
end
isi_norm = cat(1,isi_norm,B_bins);

save -append data_V1_all isi_norm bins2;


%% Bin log10 ISIs
disp('Binning log10 ISIs...')
bin_max = 3.6;
db = 0.05;
bins = 0:db:bin_max;
log_bins2 = db/2:db:bin_max-db/2;
N_bins = length(bins) - 1;

log_isi_norm = [];
load data_V1.mat
N_cells = length(isi_all);
B_bins = zeros(N_cells,N_bins);
for nc = 1:N_cells
    for nb = 1:N_bins
        id = find(log10(isi_all{nc}) > bins(nb) & log10(isi_all{nc}) <= bins(nb+1));
        B_bins(nc,nb) = length(id);
    end
    B_bins(nc,:) = B_bins(nc,:)/sum(B_bins(nc,:));
end
log_isi_norm = cat(1,log_isi_norm,B_bins);
save -append data_V1_all log_isi_norm log_bins2;

%% Pre and Post ISI
disp('Pre and Post ISI...')
isi_pre_all = [];
isi_post_all = [];
load data_V1.mat
N_cells = length(isi_all);
for nc = 1:N_cells
    isi_pre_all = cat(2,isi_pre_all,isi_pre{nc});
    isi_post_all = cat(2,isi_post_all,isi_post{nc});
end
save -append data_V1_all isi_pre_all isi_post_all;


end

function calculate_isi(N_cells,N_trials,N_stim1,N_stim2,datdir,samp_freq,refract_period,burst_thr)

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
        eval(['load ' datdir{nc} 'trial_' num2str(nt) '.mat data;'])
        
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
eval(['save data_V1 isi_all isi_pre isi_post burst_num burst2_preISI burst2_postISI burst3_preISI burst3_postISI burst_all_postISI burst_all_post2ISI'])
end

