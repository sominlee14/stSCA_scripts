% calculates spatiotemporal spike-centered average (st-SCA)
% term spatiotemporal STA is used throughout--old naming convention 
% Arguments
% file_id: patient id + recording id
% llfp_type: whitened or unwhitened
% data_site: Columbia or Harvard (for uploading the correct Utah map
% seizure_onset, seizure_end: seizure start and end time in recording
% clip_length: temporal length of clip 
% downsample_factor: 1 = no downsampling

function fx_calculate_spatiotemporal_sta(...
    file_id, llfp_type, data_site, seizure_onset, seizure_end,...
    clip_length, downsample_factor)
    
    % Load appropriate LLFP

    if strcmp(llfp_type, 'unwhitened')
        load([file_id, '_llfp.mat'], 'data', 'metadata')
    elseif strcmp(llfp_type, 'whitened')
        load([file_id, '_whitened_llfp.mat'], 'data', 'metadata')
    end

    % load spike information
    load([file_id, '_spike_detection_results.mat'], 'spike_detection_results')

    % load appropriate utah map according to data collection site
    load([data_site, '_utah_map.mat'], 'utah_map')

    sample_rate = metadata.sample_rate;
    ds_sample_rate = sample_rate / downsample_factor;
    channel_list = metadata.channel_list;
    n_channels = size(data, 1);
    n_samples_clip = ds_sample_rate * clip_length;
    
    sta_calculation_parameters = struct(...
        'file_id', file_id,...
        'seizure_onset', seizure_onset,...
        'seizure_end', seizure_end,...
        'clip_length', clip_length,...
        'sample_rate', ds_sample_rate);

    % Allocate stuff 
    weighted_channel_stas = nan(n_channels, n_samples_clip, n_channels);
    weighted_channel_stas_even = nan(n_channels, n_samples_clip, n_channels);
    weighted_channel_stas_odd = nan(n_channels, n_samples_clip, n_channels);

    channel_weights = nan(1, n_channels);
    channel_weights_even = nan(1, n_channels);
    channel_weights_odd = nan(1, n_channels);
    
    % Loop for grabbing all the st-SCAs by channel
    n = 1; % counter for keeping track of even and odd spikes for plus-minus average
    for i = 1:n_channels
        ch = channel_list(i);

        datetime('now')
        fprintf(['Running channel ', num2str(ch), '\n'])

        ch_spike_times = spike_detection_results.spike_times{i};
        ch_spike_times = ch_spike_times(...
            ch_spike_times > seizure_onset & ch_spike_times < seizure_end);

        ch_n_spikes_in_seizure = length(ch_spike_times);
        fprintf(['Spikes in channel: ', num2str(ch_n_spikes_in_seizure), '\n'])

        downsampled_spikes = nan(n_channels, n_samples_clip, ch_n_spikes_in_seizure);

        for s = 1:ch_n_spikes_in_seizure

            spike_clip = clip_signals3(data, ch_spike_times(s), sample_rate, clip_length);
            downsampled_spikes(:,:,s) = downsample_by_average(spike_clip, downsample_factor);
        
        end

        % for keeping track of whether current channel had even or odd
        % number of spikes

        if n == 1
    
            odd_spikes = downsampled_spikes(:,:,1:2:end);
            even_spikes = downsampled_spikes(:,:,2:2:end);

        elseif n == 2

            odd_spikes = downsampled_spikes(:,:,2:2:end);
            even_spikes = downsampled_spikes(:,:,1:2:end);

        end
    
        n_odd_spikes = size(odd_spikes, 3);
        n_even_spikes = size(even_spikes, 3);

        if n_odd_spikes ~= n_even_spikes

            if n == 1
                n = 2;

            elseif n == 2
                n = 1;
            end
        end
        
        weighted_channel_stas(:,:,i) = ...
            mean(downsampled_spikes,3,'omitnan') * ch_n_spikes_in_seizure;
        channel_weights(i) = ch_n_spikes_in_seizure;

        weighted_channel_stas_odd(:,:,i) = mean(odd_spikes,3,'omitnan') * n_odd_spikes;
        channel_weights_odd(i) = n_odd_spikes;

        weighted_channel_stas_even(:,:,i) = mean(even_spikes,3,'omitnan') * n_even_spikes;
        channel_weights_even(i) = n_even_spikes;

    end
    
    % free up some memory
    clearvars data

    n_rows_utah = size(utah_map,1);
    n_cols_utah = size(utah_map, 2);
    
    n_rows_field = n_rows_utah * 2 - 1;
    n_cols_field = n_cols_utah * 2 - 1;
    
    % Allocate variables for converting SCAs into grid form 
    centered_stas = nan(n_rows_field, n_cols_field, n_samples_clip, n_channels);
    centered_stas_odd = nan(n_rows_field, n_cols_field, n_samples_clip, n_channels);
    centered_stas_even = nan(n_rows_field, n_cols_field, n_samples_clip, n_channels);
    
    % Loop for converting electrodes into grid 
    for i = 1:n_channels

        ch = channel_list(i);
        fprintf(['Centering channel ', num2str(ch), '\n'])

        % Find channel's location in utah map
        [row, col] = find(utah_map == ch);
        row_diff = 10 - row;
        col_diff = 10 - col;

        channel_centered_sta = nan(n_rows_field, n_cols_field, n_samples_clip);
        channel_centered_sta_odd = nan(n_rows_field, n_cols_field, n_samples_clip);
        channel_centered_sta_even = nan(n_rows_field, n_cols_field, n_samples_clip);
        
        
        % Translate all the spikes 
        for t = 1:n_samples_clip

            for j = 1:n_channels

                j_ch = channel_list(j);
                [j_row, j_col] = find(utah_map == j_ch);
                shifted_row = j_row + row_diff;
                shifted_col = j_col + col_diff;

                channel_centered_sta(shifted_row, shifted_col, t) = ...
                    weighted_channel_stas(j, t, i);
                channel_centered_sta_odd(shifted_row, shifted_col, t) = ...
                    weighted_channel_stas_odd(j, t, i);
                channel_centered_sta_even(shifted_row, shifted_col, t) = ...
                    weighted_channel_stas_even(j, t, i);
            end

        end

        centered_stas(:,:,:,i) = channel_centered_sta;
        centered_stas_odd(:,:,:,i) = channel_centered_sta_odd;
        centered_stas_even(:,:,:,i) = channel_centered_sta_even;

    end
    
    % Find contributions for each cell in st-SCA grid 
    grid_weights = permute(sum(centered_stas,3), [1,2,4,3])*0;
    grid_weights_odd = permute(sum(centered_stas_odd,3), [1,2,4,3])*0;
    grid_weights_even = permute(sum(centered_stas_even,3), [1,2,4,3])*0;

    for i = 1:n_channels

        grid_weights(:,:,i) = grid_weights(:,:,i) + channel_weights(i);
        grid_weights_odd(:,:,i) = grid_weights_odd(:,:,i) + channel_weights_odd(i);
        grid_weights_even(:,:,i) = grid_weights_even(:,:,i) + channel_weights_even(i);

    end

    grid_weights = sum(grid_weights,3,'omitnan');
    grid_weights_odd = sum(grid_weights_odd, 3, 'omitnan');
    grid_weights_even = sum(grid_weights_even, 3, 'omitnan');

    % final st-SCA calculation 
    spatiotemporal_sta = sum(centered_stas, 4, 'omitnan')./grid_weights;
    spatiotemporal_sta_odd = sum(centered_stas_odd, 4, 'omitnan')./grid_weights_odd;
    spatiotemporal_sta_even = sum(centered_stas_even, 4, 'omitnan')./grid_weights_even;
    spatiotemporal_sta_noise_estimate = spatiotemporal_sta_odd - spatiotemporal_sta_even;

    clip_length_ms = clip_length * 1000;

    if strcmp(llfp_type, 'unwhitened')
        save_filename = [...
            file_id,...
            '_spatiotemporal_sta_',...
            num2str(clip_length_ms), 'ms_',...
            num2str(ds_sample_rate), 'hz.mat'];
        
    elseif strcmp(llfp_type, 'whitened')
        save_filename = [...
            file_id,...
            '_spatiotemporal_sta_whitened_',...
            num2str(clip_length_ms), 'ms_',...
            num2str(ds_sample_rate), 'hz.mat'];
    end

    save(save_filename,...
        'spatiotemporal_sta', 'spatiotemporal_sta_noise_estimate',...
        'grid_weights', 'sta_calculation_parameters', 'metadata')

end
