function [X, res] = getFFTintFq(T, signal, numzeros, intFq)
% Computes Fourier transform:
% Inputs T:=Sampling time, signal.
% Outputs X:=Frequency Y:=Spectral component

signal = cat(2,signal,zeros(numzeros,1)');

Fs=1/T;
L = numel(signal);
X = Fs/2*linspace(0,1,L/2)';
if sum(signal) == 0
    Y = zeros(floor(L/2),1);
else    
    Y = fft(signal ,L)/L;
    Y = Y(1:floor(L/2))';
end

res = Y(intFq);

end