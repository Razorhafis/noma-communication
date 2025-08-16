clear all;
close all;

% Parameters
%Artificial noise injection not done on this code

fs = 10e3; % Sampling frequency (10 kHz)
T = 1; % Time duration (1 second)
t = 0:1/fs:T-1/fs; % Time vector

% OFDM Parameters
N = 64; % Number of subcarriers
cp_len = 16; % Cyclic prefix length
subcarrier_spacing = fs / N; % Subcarrier spacing
f = subcarrier_spacing * (0:N-1); % Subcarrier frequencies
ofdm_symbol_duration = N / fs; % Duration of one OFDM symbol (without CP)
num_symbols = floor(T / (ofdm_symbol_duration + cp_len/fs)); % Number of OFDM symbols
samples_per_symbol = N + cp_len; % Total samples per OFDM symbol (with CP)
total_samples = num_symbols * samples_per_symbol; % Total samples
t_adjusted = 0:1/fs:(total_samples-1)/fs; % Adjusted time vector

% Generate unique sinewave signals for each user (different frequencies)
num_users = 10;
signals = zeros(total_samples, num_users); % Sinewave signals
for i = 1:num_users
    f = 500 * i; % Unique frequency for each user (0.5 kHz to 5 kHz)
    signal = sin(2 * pi * f * t_adjusted); % Sinewave
    signal_power = mean(abs(signal).^2); % Compute power of the sinewave
    signals(:, i) = signal / sqrt(signal_power); % Normalize to unit power
end

% User distances and positions
area_size = 10.0; % 10x10 meter area
uav_altitude = 10.0;
rng(12345); % Set random seed for reproducibility
uav_pos = [area_size/2, area_size/2, uav_altitude];
user_pos = zeros(num_users, 3);
dists = zeros(num_users, 1);
fprintf('User Coordinates:\n');
for i = 1:num_users
    x = rand * area_size;
    y = rand * area_size;
    user_pos(i, :) = [x, y, 0];
    dists(i) = norm(user_pos(i , :) - uav_pos);
    fprintf('User %d position: (%.5f, %.5f, %.1f), Distance: %.5f m\n', i-1, x, y, 0, dists(i));
end

% Calculate distances between users
fprintf('\nDistances Between Users:\n');
for i = 1:num_users
    for j = i+1:num_users
        dist = sqrt((user_pos(i,1) - user_pos(j,1))^2 + (user_pos(i,2) - user_pos(j,2))^2);
        fprintf('Distance between User %d and User %d: %.5f m\n', i-1, j-1, dist);
    end
end

% Path loss (simple model: 1/d^2)
path_losses = 1 ./ (dists.^2);
fprintf('\nPath Losses:\n');
for i = 1:num_users
    fprintf('User %d: Path Loss = %.5f\n', i-1, path_losses(i));
end

% Sort users by path loss (ascending) for NOMA decoding order
[~, idx_order] = sort(path_losses); % Stronger channel (higher path loss) last

% Power allocation (fixed: more power to far user)
P_total = 2; % Total power in Watts
alpha = [0.25, 0.2, 0.15, 0.12, 0.1, 0.08, 0.06, 0.05, 0.04, 0.03]; % Power allocation coefficients
powers = zeros(num_users, 1);
fprintf('\nPower Allocation:\n');
fprintf('Total Power: %.1f W\n', P_total);
for i = 1:num_users
    idx = idx_order(i);
    powers(idx) = P_total * alpha(i);
    fprintf('User %d: Power = %.3f W\n', idx-1, powers(idx));
end

% Scale signals by power
tx_signals = zeros(total_samples, num_users);
for i = 1:num_users
    tx_signals(:, i) = sqrt(powers(i)) * signals(:, i);
end

% Superimposed signal (transmitted signal)
tx_signal = sum(tx_signals, 2);

% Add noise (AWGN)
noise_power = 1e-6;
noise = sqrt(noise_power/2) * (randn(size(tx_signal)) + 1j*randn(size(tx_signal)));

% Design bandpass filters for each user
filters = cell(num_users, 1);
for i = 1:num_users
    center_freq = 500 * i; % Center frequency (0.5 kHz to 5 kHz)
    passband = 200; % Passband width (Hz)
    f_low = (center_freq - passband/2) / (fs/2); % Normalized lower frequency
    f_high = (center_freq + passband/2) / (fs/2); % Normalized upper frequency
    filters{i} = designfilt('bandpassfir', 'FilterOrder', 100, ...
        'CutoffFrequency1', f_low, 'CutoffFrequency2', f_high, ...
        'SampleRate', fs);
    fprintf('User %d Filter: Bandpass %.0f Hz to %.0f Hz\n', i-1, center_freq-passband/2, center_freq+passband/2);
end

% SINR, decoding, and data transfer
sinrs = zeros(num_users, 1);
rates = zeros(num_users, 1);
data_sent = zeros(num_users, 1);
data_received = zeros(num_users, 1);
sinr_threshold_db = 0.0; % SINR threshold for successful reception
decoded_signals = zeros(total_samples, num_users);
fprintf('\nData Transfer Details:\n');
for k = 1:num_users
    idx = idx_order(k); % User index in SIC order
    % Received signal
    rx_signal = tx_signal * path_losses(idx) + noise;
    % Apply bandpass filter to isolate user's signal
    filtered_rx_signal = filter(filters{idx}, rx_signal);
    % SIC: Decode own signal, treat later users as interference
    interference = 0;
    for i = k+1:num_users
        interference = interference + path_losses(idx) * sqrt(powers(idx_order(i))) * filter(filters{idx}, signals(:, idx_order(i)));
    end
    signal_power = mean(abs(sqrt(powers(idx)) * signals(:, idx)).^2);
    interference_power = mean(abs(interference).^2);
    sinr = (path_losses(idx)^2 * signal_power) / (interference_power + noise_power);
    sinr_db = 10 * log10(sinr);
    sinrs(idx) = sinr_db;
    rate = log2(1 + sinr);
    rates(idx) = rate;

    % Data sent (in KB)
    data_sent_bits = rate * (fs / samples_per_symbol) * T * total_samples;
    data_sent_kb = (data_sent_bits / 8) / 1024;
    data_sent(idx) = data_sent_kb;

    % Data received (based on SINR threshold)
    data_received_kb = (sinr_db >= sinr_threshold_db) * data_sent_kb;
    data_received(idx) = data_received_kb;

    % Decode signal (for plotting)
    remaining_signal = filtered_rx_signal;
    for i = 1:k-1
        remaining_signal = remaining_signal - path_losses(idx) * sqrt(powers(idx_order(i))) * filter(filters{idx}, signals(:, idx_order(i)));
    end
    decoded_signals(:, idx) = remaining_signal / (path_losses(idx) * sqrt(powers(idx)));

    fprintf('User %d: SINR = %.5f dB, Rate = %.5f bps/Hz, Data Sent = %.2f KB, Data Received = %.2f KB\n', ...
        idx-1, sinr_db, rate, data_sent_kb, data_received_kb);
end

% Summary
sum_rate = sum(rates);
total_data_sent = sum(data_sent);
total_data_received = sum(data_received);
fprintf('\nSummary:\n');
fprintf('Sum Rate: %.5f bps/Hz\n', sum_rate);
fprintf('Total Data Sent: %.2f KB\n', total_data_sent);
fprintf('Total Data Received: %.2f KB\n', total_data_received);
fprintf('Data Reception Success Rate: %.2f%%\n', (total_data_received / total_data_sent) * 100);

% Plotting Transmitted and Decoded Signals
figure;
subplot(11,1,1);
plot(t_adjusted, real(tx_signal), 'b', 'DisplayName', 'Superimposed Signal');
hold on;
plot(t_adjusted, real(tx_signals(:, 1)), 'r--', 'DisplayName', 'User 0 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 2)), 'g--', 'DisplayName', 'User 1 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 3)), 'm--', 'DisplayName', 'User 2 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 4)), 'k--', 'DisplayName', 'User 3 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 5)), 'c--', 'DisplayName', 'User 4 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 6)), 'y--', 'DisplayName', 'User 5 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 7)), 'Color', [0.5 0 0.5], 'LineStyle', '--', 'DisplayName', 'User 6 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 8)), 'Color', [0 0.5 0.5], 'LineStyle', '--', 'DisplayName', 'User 7 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 9)), 'Color', [0.5 0.5 0], 'LineStyle', '--', 'DisplayName', 'User 8 Signal (Sine)');
plot(t_adjusted, real(tx_signals(:, 10)), 'Color', [0.7 0.2 0.2], 'LineStyle', '--', 'DisplayName', 'User 9 Signal (Sine)');
title('Transmitted Superimposed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
grid on;

for u = 1:num_users
    subplot(11,1,u+1);
    plot(t_adjusted, real(decoded_signals(:, u)), 'b', 'DisplayName', sprintf('Decoded User %d', u-1));
    hold on;
    plot(t_adjusted, real(signals(:, u)), 'r--', 'DisplayName', sprintf('Original User %d', u-1));
    title(sprintf('User %d Decoded Signal', u-1));
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend;
    grid on;
end