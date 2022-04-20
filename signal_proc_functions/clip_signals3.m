function clips = clip_signals3(multichannel_signal, time_stamps, sample_rate, clip_length)

    n_samples = round(clip_length * sample_rate);
    clips = nan(size(multichannel_signal,1),  n_samples, length(time_stamps));
    i_time_stamps = round(time_stamps * sample_rate);
    
    for i = 1:length(time_stamps)
        
        i_clip_start = i_time_stamps(i) - round(n_samples/2);
        clips(:,:,i) = multichannel_signal(:, i_clip_start:(i_clip_start+n_samples-1)); 
    end
end
