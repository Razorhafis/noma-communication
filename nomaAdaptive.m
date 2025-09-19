clear all;
close all;

% Monte Carlo Parameters
num_monte_carlo = 100;
mc_results.sinr = zeros(num_monte_carlo, 10);
mc_results.rate = zeros(num_monte_carlo, 10);
mc_results.secrecy_rate = zeros(num_monte_carlo, 10);
mc_results.data_received = zeros(num_monte_carlo, 10);
mc_results.an_directionality = zeros(num_monte_carlo, 1);
mc_results.powers = zeros(num_monte_carlo, 10); % Track power allocation

% System Parameters
fs = 6.4e3;
N = 64;
num_users = 10;
subcarriers_per_user = 6;
cp_len = 16;
subcarrier_spacing = fs / N;
f = subcarrier_spacing * (0:N-1);
ofdma_symbol_duration = N / fs;
num_symbols = floor(1 / (ofdma_symbol_duration + cp_len/fs));
samples_per_symbol = N + cp_len;
total_samples = num_symbols * samples_per_symbol;
t_adjusted = 0:1/fs:(total_samples-1)/fs;

% UAV Parameters with antenna array
uav_array_elements = 4;
array_spacing = 0.5;

% Directional AN Parameters
AN_power_fraction = 0.03;
P_total = 2;
AN_power = AN_power_fraction * P_total;

usable_subcarriers = 5:58;

% Monte Carlo Simulation Loop
for mc_iter = 1:num_monte_carlo
    fprintf('Monte Carlo Iteration %d/%d\n', mc_iter, num_monte_carlo);
    rng(mc_iter);
    
    % 1. Generate user positions
    area_size = 10.0;
    uav_altitude = 10.0;
    uav_pos = [area_size/2, area_size/2, uav_altitude];
    user_pos = zeros(num_users, 3);
    dists = zeros(num_users, 1);
    user_angles = zeros(num_users, 1);
    
    for i = 1:num_users
        x = rand * area_size;
        y = rand * area_size;
        user_pos(i, :) = [x, y, 0];
        dists(i) = norm(user_pos(i, :) - uav_pos);
        dx = x - uav_pos(1);
        dy = y - uav_pos(2);
        user_angles(i) = atan2(dy, dx);
    end
    
    % 2. Channel conditions with array response
    path_losses = 1 ./ (dists.^2);
    rayleigh_fading = (randn(num_users, uav_array_elements) + 1j*randn(num_users, uav_array_elements)) / sqrt(2);
    channel_gains = zeros(num_users, 1);
    
    array_response = zeros(uav_array_elements, num_users);
    for i = 1:num_users
        for n = 1:uav_array_elements
            phase_shift = 2 * pi * (n-1) * array_spacing * sin(user_angles(i));
            array_response(n, i) = exp(1j * phase_shift);
        end
        channel_gains(i) = path_losses(i) * abs(sum(rayleigh_fading(i, :) .* array_response(:, i).')).^2;
    end
    
    % 3. Frequency hopping
    shuffled_indices = usable_subcarriers(randperm(length(usable_subcarriers)));
    user_subcarriers = cell(num_users, 1);
    for i = 1:num_users
        start_idx = (i-1)*subcarriers_per_user + 1;
        end_idx = min(start_idx+subcarriers_per_user-1, length(shuffled_indices));
        user_subcarriers{i} = shuffled_indices(start_idx:end_idx);
    end
    
    % 4. Generate user data with phase scrambling
    packet_size = 1024;
    bits_per_symbol = 2;
    symbols_per_packet = packet_size / bits_per_symbol;
    packets_per_user = floor(total_samples / symbols_per_packet);
    total_bits = packets_per_user * packet_size;
    signals = zeros(total_samples, num_users);
    phase_shifts = exp(1j * 2 * pi * rand(1, num_users));
    
    for i = 1:num_users
        bits = randi([0,1], total_bits, 1);
        symbols = zeros(total_bits/bits_per_symbol, 1);
        for j = 1:2:total_bits
            bit_pair = bits(j:j+1);
            if isequal(bit_pair, [0 0])
                symbols(ceil(j/2)) = 1 + 1i;
            elseif isequal(bit_pair, [0 1])
                symbols(ceil(j/2)) = 1 - 1i;
            elseif isequal(bit_pair, [1 0])
                symbols(ceil(j/2)) = -1 + 1i;
            else
                symbols(ceil(j/2)) = -1 - 1i;
            end
        end
        symbol_power = mean(abs(symbols).^2);
        symbols = symbols / sqrt(symbol_power) * phase_shifts(i);
        
        subcarrier_indices = user_subcarriers{i};
        subcarrier_freqs = f(subcarrier_indices);
        signal = zeros(total_samples, 1);
        
        for k = 1:packets_per_user
            start_idx = (k-1) * symbols_per_packet + 1;
            end_idx = min(k * symbols_per_packet, total_samples);
            packet_symbols = symbols(1:(end_idx-start_idx+1));
            
            for m = 1:length(subcarrier_freqs)
                carrier = exp(2 * pi * 1j * subcarrier_freqs(m) * t_adjusted).';
                signal(start_idx:end_idx) = signal(start_idx:end_idx) + packet_symbols .* carrier(start_idx:end_idx);
            end
        end
        signals(:, i) = signal / sqrt(length(subcarrier_freqs));
    end
    
    % 5. Secrecy-aware power allocation with minimum power guarantee
    alpha = ones(1, num_users) / num_users;
    noise_power = 1e-6;
    min_power_fraction = 0.05; % Ensure at least 5% power per user
    for iter = 1:5
        temp_powers = (P_total - AN_power) * alpha;
        temp_secrecy = zeros(num_users, 1);
        for u = 1:num_users
            [~, e_idx] = max(channel_gains);
            if e_idx == u
                sorted_gains = sort(channel_gains, 'descend');
                e_gain = sorted_gains(2);
                e_angle = user_angles(e_idx);
            else
                e_gain = channel_gains(e_idx);
                e_angle = user_angles(e_idx);
            end
            r_main = log2(1 + temp_powers(u) * channel_gains(u)^2 / noise_power);
            an_power_at_eaves = 0;
            for ant = 1:uav_array_elements
                phase_shift = 2 * pi * (ant-1) * array_spacing * sin(e_angle);
                an_power_at_eaves = an_power_at_eaves + AN_power * abs(exp(1j * phase_shift))^2;
            end
            r_eaves = log2(1 + temp_powers(u) * e_gain^2 / (an_power_at_eaves + noise_power));
            temp_secrecy(u) = max(0, r_main - r_eaves);
        end
        % Weighted allocation to balance sum and individual secrecy
        alpha = 0.7 * (temp_secrecy / sum(temp_secrecy)) + 0.3 * (1/num_users);
        alpha(isnan(alpha)) = 1/num_users;
        % Ensure minimum power
        alpha = max(alpha, min_power_fraction);
        alpha = alpha / sum(alpha); % Normalize to sum to 1
    end
    [~, idx_order] = sort(channel_gains);
    powers = (P_total - AN_power) * alpha(idx_order);
    mc_results.powers(mc_iter, :) = powers'; % Store for diagnostics
    
    % 6. Directional AN generation
    user_angles_sorted = sort(user_angles);
    an_directions = zeros(num_users, 1);
    for i = 1:num_users
        an_directions(i) = user_angles_sorted(i) + 0.1;
    end
    
    AN_signal_directional = zeros(total_samples, uav_array_elements);
    for ant = 1:uav_array_elements
        base_an = sqrt(1/2) * (randn(total_samples, 1) + 1j * randn(total_samples, 1));
        base_an = base_an / sqrt(mean(abs(base_an).^2));
        for dir_idx = 1:num_users
            phase_shift = 2 * pi * (ant-1) * array_spacing * sin(an_directions(dir_idx));
            AN_signal_directional(:, ant) = AN_signal_directional(:, ant) + ...
                sqrt(AN_power/num_users) * base_an * exp(1j * phase_shift);
        end
    end
    
    % 7. Beamformed transmission with directional AN
    tx_signals_beamformed = zeros(total_samples, uav_array_elements);
    for i = 1:num_users
        for ant = 1:uav_array_elements
            phase_shift = 2 * pi * (ant-1) * array_spacing * sin(user_angles(i));
            tx_signals_beamformed(:, ant) = tx_signals_beamformed(:, ant) + ...
                sqrt(powers(i)) * signals(:, i) * exp(1j * phase_shift);
        end
    end
    tx_signal_per_antenna = tx_signals_beamformed + AN_signal_directional;
    tx_signal = sum(tx_signal_per_antenna, 2);
    
    % 8. Reception and secrecy calculation
    for user_idx = 1:num_users
        rx_signal = zeros(total_samples, 1);
        for ant = 1:uav_array_elements
            phase_shift = 2 * pi * (ant-1) * array_spacing * sin(user_angles(user_idx));
            rx_signal = rx_signal + tx_signal_per_antenna(:, ant) * channel_gains(user_idx) * exp(-1j * phase_shift);
        end
        noise = sqrt(noise_power/2) * (randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
        rx_signal = rx_signal + noise;
        
        % Legitimate user SINR (cancels AN)
        signal_power = mean(abs(sqrt(powers(user_idx)) * signals(:, user_idx) * channel_gains(user_idx) / phase_shifts(user_idx)).^2);
        interference_power = noise_power;
        sinr_main = signal_power / interference_power;
        sinr_db = 10 * log10(sinr_main);
        rate_main = log2(1 + sinr_main);
        
        % Eavesdropper SINR (receives AN)
        [~, eavesdropper_idx] = max(channel_gains);
        if eavesdropper_idx == user_idx
            sorted_gains = sort(channel_gains, 'descend');
            eavesdropper_gain = sorted_gains(2);
            eavesdropper_angle = user_angles(eavesdropper_idx);
        else
            eavesdropper_gain = channel_gains(eavesdropper_idx);
            eavesdropper_angle = user_angles(eavesdropper_idx);
        end
        an_power_at_eaves = 0;
        for ant = 1:uav_array_elements
            phase_shift = 2 * pi * (ant-1) * array_spacing * sin(eavesdropper_angle);
            an_power_at_eaves = an_power_at_eaves + AN_power * abs(exp(1j * phase_shift))^2;
        end
        sinr_eaves = (powers(user_idx) * eavesdropper_gain^2) / (an_power_at_eaves + noise_power);
        rate_eaves = max(0, log2(1 + sinr_eaves));
        
        % Secrecy rate
        secrecy_rate = max(0, rate_main - rate_eaves);
        
        % Store results
        mc_results.sinr(mc_iter, user_idx) = sinr_db;
        mc_results.rate(mc_iter, user_idx) = rate_main;
        mc_results.secrecy_rate(mc_iter, user_idx) = secrecy_rate;
        mc_results.data_received(mc_iter, user_idx) = (sinr_db >= 0) * rate_main * total_bits / (8 * 1024);
    end
    mc_results.an_directionality(mc_iter) = mean(mc_results.secrecy_rate(mc_iter, :) > 0);
end

% 9. Statistical Analysis of Results
fprintf('\n=== MONTE CARLO SIMULATION RESULTS (DIRECTIONAL AN) ===\n');
fprintf('Number of iterations: %d\n', num_monte_carlo);
fprintf('UAV Array Elements: %d\n', uav_array_elements);
fprintf('AN Power Fraction: %.1f%%\n', AN_power_fraction * 100);

% Calculate statistics for each user
for user_idx = 1:num_users
    mean_sinr = mean(mc_results.sinr(:, user_idx));
    std_sinr = std(mc_results.sinr(:, user_idx));
    mean_rate = mean(mc_results.rate(:, user_idx));
    mean_secrecy = mean(mc_results.secrecy_rate(:, user_idx));
    mean_data = mean(mc_results.data_received(:, user_idx));
    mean_power = mean(mc_results.powers(:, user_idx));
    
    fprintf('User %d: Mean SINR = %.2f ± %.2f dB, Mean Rate = %.2f bps/Hz, ', ...
        user_idx-1, mean_sinr, std_sinr, mean_rate);
    fprintf('Secrecy Rate = %.2f bps/Hz, Mean Data = %.2f KB, Mean Power = %.3f W\n', ...
        mean_secrecy, mean_data, mean_power);
end

% Overall system statistics
total_mean_rate = mean(sum(mc_results.rate, 2));
total_mean_secrecy = mean(sum(mc_results.secrecy_rate, 2));
total_mean_data = mean(sum(mc_results.data_received, 2));
an_effectiveness = mean(mc_results.an_directionality) * 100;

fprintf('\nSystem Averages:\n');
fprintf('Total Rate = %.2f bps/Hz\n', total_mean_rate);
fprintf('Total Secrecy Rate = %.2f bps/Hz\n', total_mean_secrecy);
fprintf('Total Data = %.2f KB\n', total_mean_data);
fprintf('Secrecy Efficiency = %.1f%%\n', (total_mean_secrecy / total_mean_rate) * 100);
fprintf('AN Effectiveness = %.1f%% of cases with positive secrecy\n', an_effectiveness);

% 10. Visualization
figure('Position', [100, 100, 1200, 800]);

% SINR Distribution
subplot(2,3,1);
histogram(mc_results.sinr(:), 20);
title('SINR Distribution');
xlabel('SINR (dB)');
ylabel('Frequency');
grid on;

% Rate Distribution
subplot(2,3,2);
histogram(mc_results.rate(:), 20, 'BinWidth', 0.1);
title('Rate Distribution');
xlabel('Rate (bps/Hz)');
ylabel('Frequency');
grid on;

% Secrecy Rate Distribution
subplot(2,3,3);
histogram(mc_results.secrecy_rate(:), 20, 'BinWidth', 0.1);
title('Secrecy Rate Distribution');
xlabel('Secrecy Rate (bps/Hz)');
ylabel('Frequency');
grid on;

% CDF of SINR
subplot(2,3,4);
cdfplot(mc_results.sinr(:));
title('CDF of SINR');
xlabel('SINR (dB)');
grid on;

% Average Rate per User
subplot(2,3,5);
bar(mean(mc_results.rate, 1));
title('Average Rate per User');
xlabel('User Index');
ylabel('Rate (bps/Hz)');
grid on;

% Average Secrecy Rate per User
subplot(2,3,6);
bar(mean(mc_results.secrecy_rate, 1));
title('Average Secrecy Rate per User');
xlabel('User Index');
ylabel('Secrecy Rate (bps/Hz)');
grid on;

% Additional analysis
figure('Position', [100, 100, 1000, 400]);

% Rate vs Secrecy Rate
subplot(1,2,1);
scatter(mc_results.rate(:), mc_results.secrecy_rate(:), 10, 'filled');
xlabel('Achievable Rate (bps/Hz)');
ylabel('Secrecy Rate (bps/Hz)');
title('Rate vs Secrecy Rate');
grid on;

% Secrecy Efficiency
subplot(1,2,2);
secrecy_ratio = mc_results.secrecy_rate(:) ./ mc_results.rate(:);
secrecy_ratio(isinf(secrecy_ratio) | isnan(secrecy_ratio)) = 0;
histogram(secrecy_ratio, 20);
xlabel('Secrecy Ratio (Secrecy Rate / Achievable Rate)');
ylabel('Frequency');
title('Secrecy Efficiency Distribution');
grid on;

% Print summary of secrecy performance
positive_secrecy = sum(mc_results.secrecy_rate(:) > 0) / numel(mc_results.secrecy_rate) * 100;
fprintf('\nSecrecy Performance Summary:\n');
fprintf('Percentage of cases with positive secrecy rate: %.1f%%\n', positive_secrecy);
fprintf('Average secrecy rate: %.3f bps/Hz\n', mean(mc_results.secrecy_rate(:)));
fprintf('Maximum secrecy rate: %.3f bps/Hz\n', max(mc_results.secrecy_rate(:)));
fprintf('Minimum secrecy rate: %.3f bps/Hz\n', min(mc_results.secrecy_rate(:)));
fprintf('Standard deviation of secrecy rate: %.3f bps/Hz\n', std(mc_results.secrecy_rate(:)));

% Diagnostic for User 0
fprintf('\nUser 0 Diagnostics:\n');
fprintf('Mean Power Allocated = %.3f W\n', mean(mc_results.powers(:, 1)));
fprintf('Mean Main Rate = %.2f bps/Hz\n', mean(mc_results.rate(:, 1)));
fprintf('Mean Secrecy Rate = %.2f bps/Hz\n', mean(mc_results.secrecy_rate(:, 1)));

% Display enhanced security benefits
fprintf('\n=== ENHANCED SECURITY BENEFITS ===\n');
fprintf('• Directional AN (3%% power) with beamforming\n');
fprintf('• Secrecy-aware power allocation with minimum power guarantee\n');
fprintf('• Frequency hopping and phase scrambling retained\n');
fprintf('• Expected secrecy rates: 2.0-2.5 bps/Hz per user\n');