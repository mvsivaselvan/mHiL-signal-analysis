function y = prepdata(x, cutoff, fs)

% remove offset and lowpass filter the data
% x = data matrix of size # channels x # samples
% cutoff frequency for low pass filter
% fs = sampling frequency

[b, a] = butter(6, cutoff/fs, 'low');
y = filtfilt(b, a, x'); % ' since filtfilt operates on the rows
y = y';

y = y-y(:,1);
