function burst_array = classify_spikes_fn_motion_parallax(id, N_tf1, N_tf2, N_trials, N_samples)

burst_array = zeros(N_tf1, N_tf2, N_trials, N_samples);
burst_gap = 5; % ms
for ntf1 = 1:N_tf1
    for ntf2 = 1:N_tf2
        for nt = 1:N_trials
            % Go through every spike and allocate it by ISI
            if (~isempty(id{ntf1,ntf2,nt}))
                % Set spikes to 1 with ms bins
                spike_time = ceil(id{ntf1,ntf2,nt}*20);
                burst_array(ntf1,ntf2,nt,spike_time) = 1;
                % If more than 1 spike, look for bursts
                N_id = length(id{ntf1,ntf2,nt});
                isi = diff(id{ntf1,ntf2,nt});
                if (N_id >= 2)
                    burst_count = 1;
                    for nid = 1:N_id-1
                        if(isi(nid) <= burst_gap)
                            burst_count = burst_count + 1;
                            for nbc = 1:burst_count
                                spike_time = ceil(id{ntf1,ntf2,nt}(nid+2-nbc)*20);
                                burst_array(ntf1,ntf2,nt,spike_time) = burst_count;
                            end
                        else
                            burst_count = 1;
                            spike_time = ceil(id{ntf1,ntf2,nt}(nid+1)*20);
                            burst_array(ntf1,ntf2,nt,spike_time) = burst_count;
                        end
                    end
                end
            end
        end
    end
end