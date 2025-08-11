function y = prepdata(x, cutoff, fs)

% remove offset and lowpass filter the data
% x = data matrix of size # channels x # samples
% cutoff frequency for low pass filter
% fs = sampling frequency

% first remove offset, then filter to avoid end effects
y = x - x(:,1);

[b, a] = butter(6, cutoff/fs, 'low');
y = filtfilt(b, a, y'); % ' since filtfilt operates on the rows
y = y';

