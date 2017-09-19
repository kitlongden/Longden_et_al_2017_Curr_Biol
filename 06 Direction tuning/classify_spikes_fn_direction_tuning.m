function burst_array = classify_spikes_fn_direction_tuning(id, N_directions, N_trials, N_samples)

burst_array = zeros(N_directions, N_trials, N_samples);
burst_gap = 5; % ms
for nd = 1:N_directions
    for nt = 1:N_trials
        % Go through every spike and allocate it by ISI
        if (~isempty(id{nd,nt}))
            % Set spikes to 1 with ms bins
            spike_time = ceil(id{nd,nt}*20);
            burst_array(nd,nt,spike_time) = 1;
            % If more than 1 spike, look for bursts
            N_id = length(id{nd,nt});
            isi = diff(id{nd,nt});
            if (N_id >= 2)
                burst_count = 1;
                for nid = 1:N_id-1
                    if(isi(nid) <= burst_gap)
                        burst_count = burst_count + 1;
                        for nbc = 1:burst_count
                            spike_time = ceil(id{nd,nt}(nid+2-nbc)*20);
                            burst_array(nd,nt,spike_time) = burst_count;
                        end
                    else
                        burst_count = 1;
                        spike_time = ceil(id{nd,nt}(nid+1)*20);
                        burst_array(nd,nt,spike_time) = burst_count;
                    end
                end
            end
        end
    end
end