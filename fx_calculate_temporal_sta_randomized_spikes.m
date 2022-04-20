function fx_calculate_temporal_sta_randomized_spikes(index, tb, llfp_type,...
    clip_length, downsample_factor)

    file_id = tb.file_id{index}
    seizure_onset = tb.seizure_onset(index)
    seizure_end = tb.seizure_end(index)

    % load appropriate LLFP file (whitened vs. unwhitened)
    if strcmp(llfp_type, 'unwhitened')
        load([file_id, '_llfp.mat'], 'data', 'metadata')
    
    elseif strcmp(llfp_type, 'whitened')
        load([file_id, '_whitened_llfp.mat'], 'data', 'metadata')
    end
    
    % load spike times
    load([file_id, '_randomize_spike_times_results.mat'], 'randomize_spike_times_results')
    
    sample_rate = metadata.sample_rate;
    ds_sample_rate = sample_rate/downsample_factor;
    
    sta_calculation_parameters = struct(...
        'file_id', file_id,...
        'seizure_onset', seizure_onset,...
        'seizure_end', seizure_end,...
        'clip_length', clip_length,...
        'sample_rate', ds_sample_rate);
    
    n_samples_clip = clip_length * ds_sample_rate;
    n_channels = size(data, 1);
    channel_list = metadata.channel_list;

    % Allcate variables 
    weighted_channel_stas = nan(n_channels, n_samples_clip);
    weighted_channel_stas_even = nan(n_channels, n_samples_clip);
    weighted_channel_stas_odd = nan(n_channels, n_samples_clip);
    n_spikes_in_seizure = nan(n_channels, 1);
    n_spikes_in_seizure_odd = nan(n_channels, 1);
    n_spikes_in_seizure_even = nan(n_channels, 1);

    n = 1; % counter for keeping track of even/odd spikes

    data_meaned_across_channels = mean(data, 1, 'omitnan');

    for i = 1:n_channels
        
        % Get channel number
        ch = channel_list(i);
        fprintf('Running channel %d\n', round(ch))

        % find spike times that happen within seizure
        ch_spike_times = randomize_spike_times_results.randomized_spike_times{i};
        ch_spike_times = ch_spike_times(...
            ch_spike_times > seizure_onset & ch_spike_times < seizure_end);

        ch_n_spikes = randomize_spike_times_results.n_spikes(i)
        ch_n_spikes_in_seizure = length(ch_spike_times)

        downsampled_spikes = nan(ch_n_spikes_in_seizure, n_samples_clip);

        for j = 1:ch_n_spikes_in_seizure

            spike_clip = clip_signals3(...
                data_meaned_across_channels, ch_spike_times(j), sample_rate, clip_length);

            downsampled_spikes(j,:) = downsample_by_average(spike_clip, downsample_factor);

        end
        
        if n == 1
            odd_spikes = downsampled_spikes(1:2:end, :);
            even_spikes = downsampled_spikes(2:2:end, :);

        elseif n == 2
            odd_spikes = downsampled_spikes(2:2:end, :);
            even_spikes = downsampled_spikes(1:2:end, :);
        end

        n_odd_spikes = size(odd_spikes, 1);
        n_even_spikes = size(even_spikes,1);

        if n_odd_spikes ~= n_even_spikes
                
            if n == 1
                n = 2;

            elseif n == 2
                n = 1;

            end

        end
        

        weighted_channel_stas(i,:) = ...
            mean(downsampled_spikes, 1) * ch_n_spikes_in_seizure;

        weighted_channel_stas_odd(i,:) = ...
            mean(odd_spikes, 1) * n_odd_spikes;

        weighted_channel_stas_even(i,:) = ...
            mean(even_spikes, 1) * n_even_spikes;

        n_spikes_in_seizure(i) = ch_n_spikes_in_seizure;
        n_spikes_in_seizure_odd(i) = n_odd_spikes;
        n_spikes_in_seizure_even(i) = n_even_spikes;

    end

    temporal_sta = sum(weighted_channel_stas, 1, 'omitnan') / sum(n_spikes_in_seizure, 'omitnan');
    temporal_sta_odd = sum(weighted_channel_stas_odd, 1, 'omitnan') / sum(n_spikes_in_seizure_odd, 'omitnan');
    temporal_sta_even = sum(weighted_channel_stas_even, 1, 'omitnan') / sum(n_spikes_in_seizure_even, 'omitnan');

    temporal_sta_noise_estimate = temporal_sta_odd - temporal_sta_even;

    clip_length_ms = clip_length * 1000;

    if strcmp(llfp_type, 'unwhitened')
        save_filename = [...
            file_id,...
            '_temporal_sta_randomized_spikes_', ...
            num2str(clip_length_ms), 'ms_',...
            num2str(ds_sample_rate), 'hz.mat'];
    elseif strcmp(llfp_type, 'whitened')
        save_filename = [...
            file_id,...
            '_temporal_sta_whitened_randomized_spikes_', ...
            num2str(clip_length_ms), 'ms_',...
            num2str(ds_sample_rate), 'hz.mat'];
    end
    

    save(save_filename, ...
        'temporal_sta', 'temporal_sta_noise_estimate',...
        'n_spikes_in_seizure',...
        'sta_calculation_parameters', 'metadata')

end




















