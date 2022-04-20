% Downsample signal by given factor by taking average of factor # of
% neighboring samples 
%
% Last update: October 24, 2018 (SL)

function signals_downsampled = downsample_by_average(signals, downsample_factor)
    
    n_signals = size(signals,1); % Number of signals fed in argument
    n_samples = size(signals,2); % length of signal in samples
    
    % Allocate space for downsampled version of signals & reconstructions
    signals_downsampled = zeros(n_signals, ceil(n_samples/downsample_factor));
    
    for i = 1:n_signals
        
        signal_raw = signals(i,:);

        
        % Add extra NaN to end of signal to make length multiple of
        % downsample factor for reshape()
        r = rem(n_samples, downsample_factor); 
        
        if r ~= 0
            signal_raw = ...
                [signal_raw, ...
                repelem(NaN, downsample_factor - r)];
        end
        
        % Reshape signal into matrix w/ number of rows = downsample factor
        signal_reshaped = reshape(signal_raw, downsample_factor, []);
        
        % Take mean of each column resulting in downsampled signal
        signals_downsampled(i,:) = mean(signal_reshaped, 1, 'omitnan');
       
    end


end
