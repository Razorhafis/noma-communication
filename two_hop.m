clear all;
close all;

% Two-Hop NOMA System Parameters
num_far_users = 4;  % Far users (weaker channels)
num_symbols = 1000;  % Number of data symbols to transmit
bits_per_symbol = 2;  % QPSK
total_bits = num_symbols * bits_per_symbol;

% Power Parameters
P_ground = 2;  % Ground station power
P_uav = 2;     % UAV transmission power
AN_power_fraction = 0.03;
AN_power = AN_power_fraction * P_uav;
noise_power = 1e-6;

% NOMA Parameters for UAV-to-Users link
num_noma_groups = 2;
users_per_group = 2;

fprintf('=== TWO-HOP NOMA SYSTEM: Ground Station → UAV → Far Users ===\n');
fprintf('Transmitting %d QPSK symbols (%d bits) per user\n', num_symbols, total_bits);

% 1. Channel conditions
% Ground Station to UAV (UAV is near user with good channel)
gs_to_uav_gain = 0.95;  % Strong channel

% UAV to Far Users (weaker channels)
uav_to_user_gains = [0.3, 0.2, 0.4, 0.25];  % Far users have weaker channels

fprintf('\nChannel Gains:\n');
fprintf('  Ground Station → UAV: %.2f\n', gs_to_uav_gain);
fprintf('  UAV → Far Users: %s\n', mat2str(uav_to_user_gains, 2));

% 2. NOMA User Grouping for UAV-to-Users link
[~, user_ordered_indices] = sort(uav_to_user_gains, 'descend');
noma_groups = cell(num_noma_groups, 1);
for g = 1:num_noma_groups
    start_idx = (g-1)*users_per_group + 1;
    end_idx = min(g*users_per_group, num_far_users);
    noma_groups{g} = user_ordered_indices(start_idx:end_idx);
end

fprintf('\nNOMA Groups (UAV to Far Users):\n');
for g = 1:num_noma_groups
    fprintf('  Group %d: Users %s\n', g, mat2str(noma_groups{g}));
end

% 3. NOMA Power Allocation for UAV-to-Users transmission
powers = zeros(num_far_users, 1);
for g = 1:num_noma_groups
    group_users = noma_groups{g};
    if length(group_users) == 2
        strong_user = group_users(1);  % Better channel in this group
        weak_user = group_users(2);    % Worse channel in this group
        powers(weak_user) = 0.7 * (P_uav - AN_power) / num_noma_groups;
        powers(strong_user) = 0.3 * (P_uav - AN_power) / num_noma_groups;
    end
end

fprintf('\nUAV Power Allocation to Far Users:\n');
for i = 1:num_far_users
    fprintf('  Far User %d: %.4f W\n', i, powers(i));
end

% 4. GENERATE REAL USER DATA BITS at Ground Station
fprintf('\n=== GROUND STATION TRANSMISSION ===\n');
fprintf('Generating user data bits at Ground Station...\n');

user_bits = cell(num_far_users, 1);
user_symbols = cell(num_far_users, 1);

for i = 1:num_far_users
    % Generate random bits at Ground Station
    user_bits{i} = randi([0, 1], total_bits, 1);
    
    % QPSK Modulation at Ground Station
    symbols = zeros(num_symbols, 1);
    for j = 1:num_symbols
        bit1 = user_bits{i}((j-1)*2 + 1);
        bit2 = user_bits{i}((j-1)*2 + 2);
        
        if bit1 == 0 && bit2 == 0
            symbols(j) = 1 + 1j;
        elseif bit1 == 0 && bit2 == 1
            symbols(j) = 1 - 1j;
        elseif bit1 == 1 && bit2 == 0
            symbols(j) = -1 + 1j;
        else
            symbols(j) = -1 - 1j;
        end
    end
    user_symbols{i} = symbols;
end

% Display original bits from Ground Station
fprintf('\nOriginal Transmitted Bits from Ground Station (First 20 bits):\n');
for i = 1:num_far_users
    fprintf('Far User %d: %s\n', i, mat2str(user_bits{i}(1:20)'));
end

% 5. FIRST HOP: Ground Station to UAV transmission
fprintf('\n=== FIRST HOP: Ground Station → UAV ===\n');
fprintf('Ground Station transmitting to UAV...\n');

% Ground Station uses simple transmission (no NOMA to UAV)
gs_transmit_power_per_user = P_ground / num_far_users;

% UAV receives signals from Ground Station
uav_received_symbols = zeros(num_symbols, num_far_users);

for i = 1:num_far_users
    % Each user's signal transmitted separately from Ground Station
    user_signal = sqrt(gs_transmit_power_per_user) * user_symbols{i};
    
    % Add channel effects and noise
    channel_signal = user_signal * gs_to_uav_gain;
    noise = sqrt(noise_power/2) * (randn(num_symbols, 1) + 1j*randn(num_symbols, 1));
    
    uav_received_symbols(:, i) = channel_signal + noise;
end

% 6. UAV PROCESSING: Decode all user data
fprintf('\nUAV decoding received data from Ground Station...\n');
uav_decoded_bits = cell(num_far_users, 1);
uav_decoding_errors = zeros(num_far_users, 1);

for i = 1:num_far_users
    % UAV decodes each user's data (simple QPSK demodulation)
    decoded_symbols = zeros(num_symbols, 1);
    for sym_idx = 1:num_symbols
        symbol = uav_received_symbols(sym_idx, i);
        real_part = real(symbol);
        imag_part = imag(symbol);
        
        % QPSK decision at UAV
        if real_part >= 0 && imag_part >= 0
            decoded_symbols(sym_idx) = 1 + 1j;
        elseif real_part >= 0 && imag_part < 0
            decoded_symbols(sym_idx) = 1 - 1j;
        elseif real_part < 0 && imag_part >= 0
            decoded_symbols(sym_idx) = -1 + 1j;
        else
            decoded_symbols(sym_idx) = -1 - 1j;
        end
    end
    
    % Convert symbols back to bits at UAV
    user_decoded_bits = zeros(total_bits, 1);
    for sym_idx = 1:num_symbols
        symbol = decoded_symbols(sym_idx);
        
        if symbol == (1 + 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 0;
            user_decoded_bits((sym_idx-1)*2 + 2) = 0;
        elseif symbol == (1 - 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 0;
            user_decoded_bits((sym_idx-1)*2 + 2) = 1;
        elseif symbol == (-1 + 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 1;
            user_decoded_bits((sym_idx-1)*2 + 2) = 0;
        else
            user_decoded_bits((sym_idx-1)*2 + 1) = 1;
            user_decoded_bits((sym_idx-1)*2 + 2) = 1;
        end
    end
    
    uav_decoded_bits{i} = user_decoded_bits;
    uav_decoding_errors(i) = sum(user_decoded_bits ~= user_bits{i});
    
    fprintf('  Far User %d: %d decoding errors at UAV\n', i, uav_decoding_errors(i));
end

% 7. SECOND HOP: UAV to Far Users using NOMA
fprintf('\n=== SECOND HOP: UAV → Far Users (NOMA Transmission) ===\n');
fprintf('UAV performing NOMA superposition coding...\n');

% UAV re-modulates the decoded data for NOMA transmission
uav_remodulated_symbols = cell(num_far_users, 1);
for i = 1:num_far_users
    symbols = zeros(num_symbols, 1);
    bits = uav_decoded_bits{i};  % Use decoded bits (may have errors)
    
    for j = 1:num_symbols
        bit1 = bits((j-1)*2 + 1);
        bit2 = bits((j-1)*2 + 2);
        
        if bit1 == 0 && bit2 == 0
            symbols(j) = 1 + 1j;
        elseif bit1 == 0 && bit2 == 1
            symbols(j) = 1 - 1j;
        elseif bit1 == 1 && bit2 == 0
            symbols(j) = -1 + 1j;
        else
            symbols(j) = -1 - 1j;
        end
    end
    uav_remodulated_symbols{i} = symbols;
end

% NOMA Superposition Coding at UAV
superimposed_symbols = zeros(num_symbols, num_noma_groups);
for g = 1:num_noma_groups
    group_users = noma_groups{g};
    fprintf('  NOMA Group %d: ', g);
    for u = 1:length(group_users)
        user_idx = group_users(u);
        superimposed_symbols(:, g) = superimposed_symbols(:, g) + ...
            sqrt(powers(user_idx)) * uav_remodulated_symbols{user_idx};
        fprintf('User%d(%.3fW) + ', user_idx, powers(user_idx));
    end
    fprintf('\b\b  \n');
end

% Add Artificial Noise at UAV
fprintf('Adding Artificial Noise at UAV...\n');
AN_symbols = sqrt(AN_power/2) * (randn(num_symbols, 1) + 1j*randn(num_symbols, 1));

% 8. Far Users Reception with SIC
fprintf('\nFar Users receiving and decoding with SIC...\n');
far_user_received_symbols = cell(num_far_users, 1);
far_user_decoded_bits = cell(num_far_users, 1);
far_user_errors = zeros(num_far_users, 1);

for user_idx = 1:num_far_users
    fprintf('\n--- Far User %d Reception from UAV ---\n', user_idx);
    
    % Find user's NOMA group
    user_group = 0;
    for g = 1:num_noma_groups
        if ismember(user_idx, noma_groups{g})
            user_group = g;
            break;
        end
    end
    
    % Create received signal at Far User
    rx_signal = zeros(num_symbols, 1);
    
    % Desired NOMA group signal
    rx_signal = rx_signal + superimposed_symbols(:, user_group) * uav_to_user_gains(user_idx);
    
    % Interference from other NOMA groups
    for g = 1:num_noma_groups
        if g ~= user_group
            rx_signal = rx_signal + superimposed_symbols(:, g) * uav_to_user_gains(user_idx);
        end
    end
    
    % Add Artificial Noise and channel noise
    rx_signal = rx_signal + AN_symbols * uav_to_user_gains(user_idx);
    noise = sqrt(noise_power/2) * (randn(num_symbols, 1) + 1j*randn(num_symbols, 1));
    rx_signal = rx_signal + noise;
    
    far_user_received_symbols{user_idx} = rx_signal;
    
    % SIC DECODING PROCESS at Far User
    fprintf('  Performing SIC decoding...\n');
    current_signal = rx_signal;
    group_users = noma_groups{user_group};
    
    % Sort users in this group by channel strength
    [~, sic_order] = sort(uav_to_user_gains(group_users), 'descend');
    
    fprintf('  SIC order in group: ');
    for k = 1:length(sic_order)
        fprintf('User%d ', group_users(sic_order(k)));
    end
    fprintf('\n');
    
    % SIC: Decode and subtract stronger users
    decoded_symbols = zeros(num_symbols, 1);
    
    for k = 1:length(sic_order)
        current_user = group_users(sic_order(k));
        
        if current_user == user_idx
            % Decode own signal
            fprintf('  Decoding own signal...\n');
            
            % QPSK demodulation
            for sym_idx = 1:num_symbols
                symbol = current_signal(sym_idx);
                real_part = real(symbol);
                imag_part = imag(symbol);
                
                if real_part >= 0 && imag_part >= 0
                    decoded_symbols(sym_idx) = 1 + 1j;
                elseif real_part >= 0 && imag_part < 0
                    decoded_symbols(sym_idx) = 1 - 1j;
                elseif real_part < 0 && imag_part >= 0
                    decoded_symbols(sym_idx) = -1 + 1j;
                else
                    decoded_symbols(sym_idx) = -1 - 1j;
                end
            end
            break;
            
        else
            % Decode and subtract stronger user
            fprintf('  Subtracting User%d interference...\n', current_user);
            
            % Decode stronger user's symbols
            temp_decoded = zeros(num_symbols, 1);
            for sym_idx = 1:num_symbols
                symbol = current_signal(sym_idx);
                real_part = real(symbol);
                imag_part = imag(symbol);
                
                if real_part >= 0 && imag_part >= 0
                    temp_decoded(sym_idx) = 1 + 1j;
                elseif real_part >= 0 && imag_part < 0
                    temp_decoded(sym_idx) = 1 - 1j;
                elseif real_part < 0 && imag_part >= 0
                    temp_decoded(sym_idx) = -1 + 1j;
                else
                    temp_decoded(sym_idx) = -1 - 1j;
                end
            end
            
            % Subtract reconstructed signal
            reconstructed = sqrt(powers(current_user)) * temp_decoded * uav_to_user_gains(user_idx);
            current_signal = current_signal - reconstructed;
        end
    end
    
    % Convert symbols back to bits
    user_decoded_bits = zeros(total_bits, 1);
    for sym_idx = 1:num_symbols
        symbol = decoded_symbols(sym_idx);
        
        if symbol == (1 + 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 0;
            user_decoded_bits((sym_idx-1)*2 + 2) = 0;
        elseif symbol == (1 - 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 0;
            user_decoded_bits((sym_idx-1)*2 + 2) = 1;
        elseif symbol == (-1 + 1j)
            user_decoded_bits((sym_idx-1)*2 + 1) = 1;
            user_decoded_bits((sym_idx-1)*2 + 2) = 0;
        else
            user_decoded_bits((sym_idx-1)*2 + 1) = 1;
            user_decoded_bits((sym_idx-1)*2 + 2) = 1;
        end
    end
    
    far_user_decoded_bits{user_idx} = user_decoded_bits;
    far_user_errors(user_idx) = sum(user_decoded_bits ~= user_bits{user_idx});
    
    fprintf('  Final: %d errors, BER: %.6f\n', far_user_errors(user_idx), far_user_errors(user_idx)/total_bits);
end

% 9. EAVESDROPPER ANALYSIS (trying to intercept UAV-to-Users transmission)
fprintf('\n=== EAVESDROPPER ANALYSIS ===\n');
eaves_gain = 0.35;  % Eavesdropper has weak channel to UAV

% Eavesdropper receives UAV transmission
eaves_signal = zeros(num_symbols, 1);
for g = 1:num_noma_groups
    eaves_signal = eaves_signal + superimposed_symbols(:, g) * eaves_gain;
end
eaves_signal = eaves_signal + AN_symbols * eaves_gain;
eaves_noise = sqrt(noise_power/2) * (randn(num_symbols, 1) + 1j*randn(num_symbols, 1));
eaves_signal = eaves_signal + eaves_noise;

% Eavesdropper tries to decode (no SIC capability)
fprintf('Eavesdropper decoding attempts (no SIC):\n');
eaves_bers = zeros(num_far_users, 1);

for user_idx = 1:num_far_users
    eaves_decoded_bits = zeros(total_bits, 1);
    for sym_idx = 1:num_symbols
        symbol = eaves_signal(sym_idx);
        real_part = real(symbol);
        imag_part = imag(symbol);
        
        if real_part >= 0 && imag_part >= 0
            eaves_decoded_bits((sym_idx-1)*2 + 1) = 0;
            eaves_decoded_bits((sym_idx-1)*2 + 2) = 0;
        elseif real_part >= 0 && imag_part < 0
            eaves_decoded_bits((sym_idx-1)*2 + 1) = 0;
            eaves_decoded_bits((sym_idx-1)*2 + 2) = 1;
        elseif real_part < 0 && imag_part >= 0
            eaves_decoded_bits((sym_idx-1)*2 + 1) = 1;
            eaves_decoded_bits((sym_idx-1)*2 + 2) = 0;
        else
            eaves_decoded_bits((sym_idx-1)*2 + 1) = 1;
            eaves_decoded_bits((sym_idx-1)*2 + 2) = 1;
        end
    end
    
    eaves_errors = sum(eaves_decoded_bits ~= user_bits{user_idx});
    eaves_ber = eaves_errors / total_bits;
    eaves_bers(user_idx) = eaves_ber;
    
    fprintf('  Far User %d: %d errors, BER: %.6f\n', user_idx, eaves_errors, eaves_ber);
end

% 10. COMPREHENSIVE RESULTS
fprintf('\n=== TWO-HOP NOMA SYSTEM RESULTS ===\n');

fprintf('\nFirst Hop (Ground Station → UAV):\n');
fprintf('User | UAV Decoding Errors | UAV BER\n');
fprintf('-----|---------------------|--------\n');
for i = 1:num_far_users
    fprintf('  %d  |         %4d         | %.6f\n', i, uav_decoding_errors(i), uav_decoding_errors(i)/total_bits);
end

fprintf('\nSecond Hop (UAV → Far Users):\n');
fprintf('User | Final Errors | Final BER | Eaves BER | Security Advantage\n');
fprintf('-----|--------------|-----------|-----------|-------------------\n');
for i = 1:num_far_users
    final_ber = far_user_errors(i)/total_bits;
    eaves_ber = eaves_bers(i);
    if final_ber > 0
        security_advantage = eaves_ber / final_ber;
    else
        security_advantage = inf;
    end
    fprintf('  %d  |     %4d     | %8.6f | %8.6f |       %6.2fx\n', ...
        i, far_user_errors(i), final_ber, eaves_ber, security_advantage);
end

% Visualization
figure('Position', [100, 100, 1400, 900]);

% System architecture
subplot(3,4,1);
plot(0, 0, 'ks', 'MarkerSize', 15, 'LineWidth', 3); hold on; % Ground Station
plot(1, 1, 'r^', 'MarkerSize', 15, 'LineWidth', 3); % UAV
for i = 1:num_far_users
    plot(2, i*0.5, 'bo', 'MarkerSize', 10, 'LineWidth', 2); % Far Users
    text(2.1, i*0.5, sprintf('User %d', i));
end
plot(3, 2.5, 'mx', 'MarkerSize', 15, 'LineWidth', 3); % Eavesdropper
text(3.1, 2.5, 'Eavesdropper');
plot([0 1], [0 1], 'k-', 'LineWidth', 2);
plot([1 2], [1 0.5], 'b-', 'LineWidth', 1);
plot([1 2], [1 1.0], 'b-', 'LineWidth', 1);
plot([1 2], [1 1.5], 'b-', 'LineWidth', 1);
plot([1 2], [1 2.0], 'b-', 'LineWidth', 1);
plot([1 3], [1 2.5], 'm--', 'LineWidth', 2);
text(0, -0.2, 'Ground Station', 'HorizontalAlignment', 'center');
text(1, 1.2, 'UAV (Relay)', 'HorizontalAlignment', 'center');
axis([-0.5 3.5 -0.5 3]);
title('Two-Hop NOMA System Architecture');
grid on;

% BER comparison
subplot(3,4,2);
uav_bers = uav_decoding_errors / total_bits;
final_bers = far_user_errors / total_bits;
bar([uav_bers; final_bers; eaves_bers]');
legend('UAV BER', 'Final BER', 'Eaves BER', 'Location', 'northeast');
title('Bit Error Rate Comparison');
xlabel('Far User Index');
ylabel('BER');
grid on;

% Channel gains
subplot(3,4,3);
gains = [gs_to_uav_gain, uav_to_user_gains, eaves_gain];
bar(gains);
title('Channel Gains');
xlabel('Link (1:GS→UAV, 2-5:UAV→Users, 6:Eaves)');
ylabel('Gain');
grid on;

% Power allocation
subplot(3,4,4);
bar(powers);
title('UAV Power Allocation to Far Users');
xlabel('Far User Index');
ylabel('Power (W)');
grid on;

% Bit comparisons for selected users
for i = 1:3
    subplot(3,4,4+i);
    bit_range = 1:min(30, total_bits);
    plot(bit_range, user_bits{i}(bit_range), 'bo-', 'LineWidth', 2, 'MarkerSize', 4);
    hold on;
    plot(bit_range, far_user_decoded_bits{i}(bit_range), 'rx--', 'LineWidth', 1, 'MarkerSize', 4);
    title(sprintf('Far User %d: Original vs Final', i));
    xlabel('Bit Index');
    ylabel('Bit Value');
    legend('Original', 'Final');
    grid on;
end

% Error progression
subplot(3,4,8);
errors = [uav_decoding_errors, far_user_errors]';
bar(errors);
title('Error Progression Through Hops');
xlabel('Far User Index');
ylabel('Number of Errors');
legend('UAV Errors', 'Final Errors', 'Location', 'northeast');
grid on;

% Security advantage
subplot(3,4,9);
security_ratios = zeros(num_far_users, 1);
for i = 1:num_far_users
    if final_bers(i) > 0
        security_ratios(i) = eaves_bers(i) / final_bers(i);
    else
        security_ratios(i) = 50; % Cap for display
    end
end
bar(security_ratios);
title('Security Advantage Ratio');
xlabel('Far User Index');
ylabel('Eaves BER / Final BER');
grid on;

% NOMA groups
subplot(3,4,10);
group_visual = zeros(num_far_users, 1);
for g = 1:num_noma_groups
    group_users = noma_groups{g};
    for u = 1:length(group_users)
        group_visual(group_users(u)) = g;
    end
end
bar(group_visual);
title('NOMA Group Assignment');
xlabel('Far User Index');
ylabel('Group Number');
grid on;

% Constellation at UAV
subplot(3,4,11);
plot(real(uav_received_symbols(:,1)), imag(uav_received_symbols(:,1)), '.', 'MarkerSize', 1);
title('UAV: Received Constellation (User 1)');
xlabel('In-phase');
ylabel('Quadrature');
grid on;
axis equal;

% Constellation at Far User
subplot(3,4,12);
plot(real(far_user_received_symbols{1}), imag(far_user_received_symbols{1}), '.', 'MarkerSize', 1);
title('Far User 1: Received Constellation');
xlabel('In-phase');
ylabel('Quadrature');
grid on;
axis equal;

fprintf('\n=== SYSTEM SUMMARY ===\n');
fprintf('Two-Hop Architecture: Ground Station → UAV → Far Users\n');
fprintf('- UAV acts as relay with good channel from Ground Station\n');
fprintf('- UAV uses NOMA to serve multiple far users simultaneously\n');
fprintf('- Far users perform SIC to decode their signals\n');
fprintf('- Eavesdropper cannot intercept due to NOMA interference + AN\n');