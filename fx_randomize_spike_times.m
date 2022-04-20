% function for randomizing spike times 
% Spike times are randomized within and across MEA channels 
% Arguments: 
% tb = table containing all the metadata
% index = row of tb to run (i.e. recording number)

function fx_randomize_spike_times(tb, index)
    
    % Get all the relevant metadata for recording 
    file_id = tb.file_id{index}
    seizure_onset = tb.seizure_onset(index)
    seizure_end = tb.seizure_end(index)
    seizure_length = seizure_end - seizure_onset
    total_n_spikes_in_seizure = tb.total_n_spikes_in_seizure(index)

    load([file_id, '_spike_detection_results.mat'], ...
    'metadata', 'spike_detection_results')

    sample_rate = metadata.sample_rate;
    n_channels = height(spike_detection_results);

    n_time_samples_seizure = seizure_length * sample_rate;
    i_randomized_spikes = randperm(n_channels * n_time_samples_seizure, total_n_spikes_in_seizure);

    
    % Initialize randomized raster 
    randomized_raster = zeros(n_channels, n_time_samples_seizure);

    % random indices for spikes
    i_randomized_spikes = randperm(n_channels * n_time_samples_seizure, total_n_spikes_in_seizure);
    
    % set random indices to 1
    randomized_raster(i_randomized_spikes) = 1;

    randomized_spike_times = cell(n_channels,1);
    n_spikes = nan(n_channels,1);

    % convert indices to times, store in cell structure
    for i = 1:n_channels
        
        % get spike times for a channel
        randomized_spike_times{i} = (find(randomized_raster(i,:))/sample_rate)+seizure_onset;
        n_spikes(i) = sum(randomized_raster(i,:));

    end
    
    randomize_spike_times_results = table(randomized_spike_times, n_spikes);


    fprintf(['n_spikes in original spike detection: ', ...
        num2str(total_n_spikes_in_seizure), '\n'])
    fprintf(['n_spikes in randomized raster: ',...
        num2str(sum(n_spikes, 'omitnan')), '\n'])

    save_filename = [file_id, '_randomize_spike_times_results.mat'];
    save(save_filename, 'randomize_spike_times_results',...
        'spike_detection_results', 'metadata', '-v7.3')



end
