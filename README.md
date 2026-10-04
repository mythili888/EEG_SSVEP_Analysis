# EEG_SSVEP_Analysis

% EEG SSVEP Analysis %

close all; clear; clc;

data_root = '/MATLAB Drive/EEG Data';

% The sampling rate as defined in the EEG Data Specification
fs = 128; % total samples/ duration : 3840/30

%% Load Data
% Loading all 5 trials for each frequency (7, 8, 9, 10 Hz)

for t = 1:5
    tmp = load(fullfile(data_root, '7Hz', ['7_' num2str(t) '.mat']));
    data_7(:,t) = double(tmp.data);

    tmp = load(fullfile(data_root, '8Hz', ['8_' num2str(t) '.mat']));
    data_8(:,t) = double(tmp.data);

    tmp = load(fullfile(data_root, '9Hz', ['9_' num2str(t) '.mat']));
    data_9(:,t) = double(tmp.data);

    tmp = load(fullfile(data_root, '10Hz', ['10_' num2str(t) '.mat']));
    data_10(:,t) = double(tmp.data);
end

disp('All data loaded.')

% Time axis - length determined from loaded data
time = (1:length(data_7(:,1)))/fs;
Signallength_time = length(data_7(:,1))/fs; %dividing by fs converts sample numbers into seconds

%% Plot Raw Signals
% Trial 1 shown as representative example to demonstrate DC offset

figure
tiledlayout(2,2)

nexttile
plot(time, data_7(:,1))
xlim([0 Signallength_time])
xlabel('Time in seconds')
ylabel('\muV')
title('7 Hz — Unfiltered EEG Signal')
grid minor

nexttile
plot(time, data_8(:,1))
xlim([0 Signallength_time])
xlabel('Time in seconds')
ylabel('\muV')
title('8 Hz — Unfiltered EEG Signal')
grid minor

nexttile
plot(time, data_9(:,1))
xlim([0 Signallength_time])
xlabel('Time in seconds')
ylabel('\muV')
title('9 Hz — Unfiltered EEG Signal')
grid minor

nexttile
plot(time, data_10(:,1))
xlim([0 Signallength_time])
xlabel('Time in seconds')
ylabel('\muV')
title('10 Hz — Unfiltered EEG Signal')
grid minor

%% Filtering %
% Trial 1 filtered as representative example to demonstrate
% effect of bandpass filtering on the raw EEG signal

figure
tiledlayout(2,2)

% [ 7 Hz ]
% Low pass filter design
filterOrder = 5;        % Order of filter
cutOffFreq = 8;         % Cutoff frequency - 1 Hz above stimulus
[b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
% Apply low pass filter
EEG_lowpass_7 = filter(b, a, data_7(:,1));

% High pass filter design
filterOrder = 4;        % Order of filter
cutOffFreq = 6;         % Cutoff frequency - 1 Hz below stimulus
[bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
% Apply high pass filter
EEG_filtered_7 = filter(bh, ah, EEG_lowpass_7);

nexttile
plot(time, EEG_filtered_7)
xlim([0 Signallength_time])
ylim([-10 10])
xlabel('Time in seconds')
ylabel('\muV')
title('7 Hz — Filtered EEG')
grid minor

% [ 8 Hz ]
filterOrder = 5;
cutOffFreq = 9;
[b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
EEG_lowpass_8 = filter(b, a, data_8(:,1));

filterOrder = 4;
cutOffFreq = 7;
[bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
EEG_filtered_8 = filter(bh, ah, EEG_lowpass_8);

nexttile
plot(time, EEG_filtered_8)
xlim([0 Signallength_time])
ylim([-10 10])
xlabel('Time in seconds')
ylabel('\muV')
title('8 Hz — Filtered EEG')
grid minor

% [ 9 Hz ]
filterOrder = 5;
cutOffFreq = 10;
[b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
EEG_lowpass_9 = filter(b, a, data_9(:,1));

filterOrder = 4;
cutOffFreq = 8;
[bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
EEG_filtered_9 = filter(bh, ah, EEG_lowpass_9);

nexttile
plot(time, EEG_filtered_9)
xlim([0 Signallength_time])
ylim([-10 10])
xlabel('Time in seconds')
ylabel('\muV')
title('9 Hz — Filtered EEG')
grid minor

% [ 10 Hz ]
filterOrder = 5;
cutOffFreq = 11;
[b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
EEG_lowpass_10 = filter(b, a, data_10(:,1));

filterOrder = 4;
cutOffFreq = 9;
[bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
EEG_filtered_10 = filter(bh, ah, EEG_lowpass_10);

nexttile
plot(time, EEG_filtered_10)
xlim([0 Signallength_time])
ylim([-10 10])
xlabel('Time in seconds')
ylabel('\muV')
title('10 Hz — Filtered EEG')
grid minor

%% FFT Computation
% Mean FFT computed across all 5 trials per frequency (7,8,9 and 10Hz)
% Peak amplitude extracted from each trial for box plot comparison

N = length(data_7(:,1)); %number of samples per trial : 30x128
f_axis = (0:N-1)*(fs/N); %FFT equation
f_axis_half = f_axis(1:floor(N/2)+1);

peak_amplitudes = zeros(4,5);

figure
tiledlayout(2,2)

% [ 7 Hz FFT ]
mean_amp_7 = zeros(length(f_axis_half), 1);
for t = 1:5
    sig = data_7(:,t);
    filterOrder = 5; cutOffFreq = 8;
    [b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
    sig_lp = filter(b, a, sig);
    filterOrder = 4; cutOffFreq = 6;
    [bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
    sig_f = filter(bh, ah, sig_lp);

    Y = fft(sig_f);
    amp = abs(Y)/N;
    amp_single = 2*amp(1:floor(N/2)+1);
    amp_single(1) = amp_single(1)/2; % DC at bin 0 should not be doubled, corrected
    mean_amp_7 = mean_amp_7 + amp_single(:);
    idx = find(f_axis_half >= 6.5 & f_axis_half <= 7.5); % 0.5 window accounts for biological variability
    peak_amplitudes(1,t) = max(amp_single(idx));
end
mean_amp_7 = mean_amp_7 / 5; %average spectrum

nexttile
plot(f_axis_half, mean_amp_7)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Amplitude (\muV)')
title('7 Hz — Mean FFT Spectrum (n=5)')
xline(7, 'r--', '7 Hz')
grid minor

% [ 8 Hz FFT ]
mean_amp_8 = zeros(length(f_axis_half), 1);
for t = 1:5
    sig = data_8(:,t);
    filterOrder = 5; cutOffFreq = 9;
    [b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
    sig_lp = filter(b, a, sig);
    filterOrder = 4; cutOffFreq = 7;
    [bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
    sig_f = filter(bh, ah, sig_lp);

    Y = fft(sig_f);
    amp = abs(Y)/N;
    amp_single = 2*amp(1:floor(N/2)+1);
    amp_single(1) = amp_single(1)/2;
    mean_amp_8 = mean_amp_8 + amp_single(:);
    idx = find(f_axis_half >= 7.5 & f_axis_half <= 8.5);
    peak_amplitudes(2,t) = max(amp_single(idx));
end
mean_amp_8 = mean_amp_8 / 5;

nexttile
plot(f_axis_half, mean_amp_8)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Amplitude (\muV)')
title('8 Hz — Mean FFT Spectrum (n=5)')
xline(8, 'r--', '8 Hz')
grid minor

% [ 9 Hz FFT ]
mean_amp_9 = zeros(length(f_axis_half), 1);
for t = 1:5
    sig = data_9(:,t);
    filterOrder = 5; cutOffFreq = 10;
    [b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
    sig_lp = filter(b, a, sig);
    filterOrder = 4; cutOffFreq = 8;
    [bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
    sig_f = filter(bh, ah, sig_lp);

    Y = fft(sig_f);
    amp = abs(Y)/N;
    amp_single = 2*amp(1:floor(N/2)+1);
    amp_single(1) = amp_single(1)/2;
    mean_amp_9 = mean_amp_9 + amp_single(:);
    idx = find(f_axis_half >= 8.5 & f_axis_half <= 9.5);
    peak_amplitudes(3,t) = max(amp_single(idx));
end
mean_amp_9 = mean_amp_9 / 5;

nexttile
plot(f_axis_half, mean_amp_9)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Amplitude (\muV)')
title('9 Hz — Mean FFT Spectrum (n=5)')
xline(9, 'r--', '9 Hz')
grid minor

% [ 10 Hz FFT ]
mean_amp_10 = zeros(length(f_axis_half), 1);
for t = 1:5
    sig = data_10(:,t);
    filterOrder = 5; cutOffFreq = 11;
    [b, a] = butter(filterOrder, cutOffFreq/(fs/2), 'low');
    sig_lp = filter(b, a, sig);
    filterOrder = 4; cutOffFreq = 9;
    [bh, ah] = butter(filterOrder, cutOffFreq/(fs/2), 'high');
    sig_f = filter(bh, ah, sig_lp);

    Y = fft(sig_f);
    amp = abs(Y)/N;
    amp_single = 2*amp(1:floor(N/2)+1);
    amp_single(1) = amp_single(1)/2;
    mean_amp_10 = mean_amp_10 + amp_single(:);
    idx = find(f_axis_half >= 9.5 & f_axis_half <= 10.5);
    peak_amplitudes(4,t) = max(amp_single(idx));
end
mean_amp_10 = mean_amp_10 / 5;

nexttile
plot(f_axis_half, mean_amp_10)
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Amplitude (\muV)')
title('10 Hz — Mean FFT Spectrum (n=5)')
xline(10, 'r--', '10 Hz')
grid minor

%% Box Plot Comparison
% Comparing peak SSVEP amplitude across all four stimulus frequencies

figure
boxplot(peak_amplitudes', ...
    'Labels', {'7 Hz','8 Hz','9 Hz','10 Hz'}, ...
    'Widths', 0.5)
xlabel('Stimulus Frequency (Hz)')
ylabel('Peak Amplitude (\muV)')
title('SSVEP Peak Amplitude Comparison Across Stimulus Frequencies')
grid minor


%% Standard Deviation Analysis for comparison between frequencies
% checks inter-trial variability, low SD indicates consistent SSVEP
% response across trials
freqs = [7, 8, 9, 10];

fprintf('\n--- STANDARD DEVIATION ---\n')
fprintf('  7 Hz:  %.4f uV\n', std(peak_amplitudes(1,:)))
fprintf('  8 Hz:  %.4f uV\n', std(peak_amplitudes(2,:)))
fprintf('  9 Hz:  %.4f uV\n', std(peak_amplitudes(3,:)))
fprintf(' 10 Hz:  %.4f uV\n', std(peak_amplitudes(4,:)))

fprintf('\n--- MEAN +/- SD ---\n')
for i = 1:4
    fprintf('  %2d Hz:  %.4f +/- %.4f uV\n', freqs(i), mean(peak_amplitudes(i,:)), std(peak_amplitudes(i,:)))
end
%% Results Summary - displays data

[max_amp, max_idx] = max(mean(peak_amplitudes, 2));
[min_amp, min_idx] = min(mean(peak_amplitudes, 2));

stim_freqs = [7, 8, 9, 10];
fprintf('\n--- RESULTS SUMMARY ---\n')
fprintf('Highest mean amplitude: %d Hz  (%.4f uV)\n', stim_freqs(max_idx), max_amp)
fprintf('Lowest  mean amplitude: %d Hz  (%.4f uV)\n', stim_freqs(min_idx), min_amp)
fprintf('\nMean peak amplitudes:\n')
fprintf('  7 Hz:  %.4f uV\n', mean(peak_amplitudes(1,:)))
fprintf('  8 Hz:  %.4f uV\n', mean(peak_amplitudes(2,:)))
fprintf('  9 Hz:  %.4f uV\n', mean(peak_amplitudes(3,:)))
fprintf(' 10 Hz:  %.4f uV\n', mean(peak_amplitudes(4,:)))
