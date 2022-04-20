function [whitened_signals, C, W] = whiten_signals(signals)

    C = cov(signals');
    [V,D] = eig(C);
    W = V * (sqrt(D)\V');
    whitened_signals = W * signals;
    whitened_signals = real(whitened_signals);





end
