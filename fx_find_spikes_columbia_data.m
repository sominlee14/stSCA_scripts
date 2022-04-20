% function for finding and saving spike times for columbia data set
% Arguments: file_id = patientID + recording number

function fx_find_spikes_columbia_data(file_id)
    
    % Load corresponding mua file 
    mua_file = [file_id, '_mua.mat'];
    load(mua_file, 'data', 'metadata')

    % run status update
    datetime('now')
    fprintf(['now running ', file_id, '\n'])
    fprintf([mua_file, ' loaded\n'])

    sample_rate = metadata.sample_rate;
    n_channels = size(data,1);
    spike_times = cell(n_channels,1);
    n_spikes = nan(n_channels,1);

    for i = 1:n_channels
        
        % for tracking run
        if mod(i,5)==0
            fprintf([num2str(i), ' of ', num2str(n_channels), ' channels done\n'])
        end
        
        % perform spike detection on single channel of data
        signal = data(i,:); 
        threshold = 5.93 * median(abs(signal(1:(15*sample_rate))));
        [~,spike_times{i}] = findpeaks(-1*signal,...
            sample_rate,...
            'MinPeakHeight', threshold,...
            'MinPeakDistance', 0.001);

        n_spikes(i) = length(spike_times{i});

    end

    spike_detection_results = table(spike_times, n_spikes);
    n_total_spikes = sum(n_spikes);
    
    % save spike detection results
    save_filename = [file_id, '_spike_detection_results.mat'];
    save(save_filename, 'spike_detection_results', 'n_total_spikes', ...
        'metadata', '-v7.3')

end
