function [X Y] = getFourier(time, signal, numOfPoints)
T = time(3)-time(2);

if numOfPoints-numel(time) > 0
        
    numOfZeros = floor((numOfPoints-numel(time))/2);
    newTime = time(numel(time)):T:( T*numOfZeros+time(numel(time)) );
    time = cat(2,time, newTime);
    newTime = (-T*numOfZeros- time(1)):T:(-time(1));
    time = cat(2,newTime,time );
    
    if numel(time)>numOfPoints
        time = time(1:numOfPoints);
    end
    
else

    disp('WARNING: Possible insufficient number of Points');

end
Fs=1/T;
L = numel(time);

% NFFT = 2^nextpow2(L);
% NFFT = L;
X = Fs/2*linspace(0,1,L/2)';
Y = fft(signal ,L)/L;
Y = Y(1:(L/2));