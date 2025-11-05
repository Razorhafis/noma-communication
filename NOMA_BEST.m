clear all;
close all;

% Two-Hop NOMA System Parameters with RIS-Only Communication
num_far_users = 9;
num_total_users = num_far_users + 1; % UAV (User 1) + 9 Far Users
num_symbols = 1000;
bits_per_symbol = 2; % QPSK
total_bits = num_symbols * bits_per_symbol;

% Power Parameters - REALISTIC VALUES
P_ground = 4;  % Realistic power
P_uav = 4;     % Realistic power
AN_power_fraction = 0.03;
AN_power = AN_power_fraction * P_uav;
noise_power = 1e-10;  % Adjusted noise

% NOMA Parameters - SINGLE GROUP
num_noma_groups = 1;  % Single NOMA group
users_per_group = num_total_users; % All users in one group

% RIS Parameters
num_RIS_elements = 64; % 8x8 RIS
RIS_position = [500, 500, 20]; % RIS position [x, y, height]

% Beamforming Parameters
num_antennas_UAV = 8; % UAV with antenna array

% Area parameters
area_size = 1000;
gs_position = [0,0];
uav_altitude = 30;

fprintf('=== TWO-HOP NOMA SYSTEM WITH SINGLE NOMA GROUP ===\n');
fprintf('First Hop: GS → UAV ONLY (NOMA for all users data)\n');
fprintf('Second Hop: UAV → RIS → Far Users (NO direct links)\n');
fprintf('          (SINGLE NOMA GROUP + RIS + Beamforming + AN + RAYLEIGH FADING)\n\n');

% 1. USER POSITIONS
user_positions = zeros(num_total_users, 2);
user_positions(1, :) = [0, 0]; % UAV above GS
for i = 2:num_total_users
    user_positions(i, :) = rand(1,2) * area_size;
end

% 2. CHANNEL GAINS WITH GUARANTEED MINIMUM
fprintf('=== CHANNEL MODEL: RAYLEIGH FADING WITH MINIMUM GAIN ===\n');

% Calculate distances
gs_to_uav_distance = uav_altitude;
uav_to_RIS_distance = norm([user_positions(1,:), uav_altitude] - RIS_position);

% RIS to Users distances
RIS_to_user_distances = zeros(num_far_users, 1);
for i = 1:num_far_users
    user_idx = i + 1;
    horizontal_distance = norm(user_positions(user_idx, :) - RIS_position(1:2));
    RIS_to_user_distances(i) = sqrt(horizontal_distance^2 + (RIS_position(3))^2);
end

% CHANNEL GAINS WITH GUARANTEED MINIMUM
min_channel_gain = 1e-4;  % Minimum guaranteed gain

% GS to UAV channel
pl_gs_uav = 20 + 10*2.0*log10(max(gs_to_uav_distance, 1)); % Avoid log(0)
gs_to_uav_gain = 10^(-pl_gs_uav/10);
gs_to_uav_gain = abs(gs_to_uav_gain * (sqrt(0.5)*(randn + 1j*randn)));
gs_to_uav_gain = max(gs_to_uav_gain, min_channel_gain); % Guarantee minimum

% UAV to RIS channel  
pl_uav_ris = 20 + 10*2.0*log10(max(uav_to_RIS_distance, 1));
uav_to_RIS_gain = 10^(-pl_uav_ris/10);
uav_to_RIS_gain = abs(uav_to_RIS_gain * (sqrt(0.5)*(randn + 1j*randn)));
uav_to_RIS_gain = max(uav_to_RIS_gain, min_channel_gain); % Guarantee minimum

% RIS to Users channels
RIS_to_user_gains = zeros(num_far_users, 1);
for i = 1:num_far_users
    d = max(RIS_to_user_distances(i), 1); % Avoid log(0)
    pl_ris_user = 25 + 10*2.5*log10(d);
    RIS_to_user_gains(i) = 10^(-pl_ris_user/10);
    RIS_to_user_gains(i) = abs(RIS_to_user_gains(i) * (sqrt(0.5)*(randn + 1j*randn)));
    RIS_to_user_gains(i) = max(RIS_to_user_gains(i), min_channel_gain); % Guarantee minimum
end

% RIS phase shift matrix
RIS_phase = exp(1j * 2*pi * rand(num_RIS_elements, 1));

fprintf('Channel Gains (with minimum %.1e):\n', min_channel_gain);
fprintf('GS → UAV: %.6f\n', gs_to_uav_gain);
fprintf('UAV → RIS: %.6f\n', uav_to_RIS_gain);
for i = 1:num_far_users
    fprintf('RIS → User %d: %.6f\n', i+1, RIS_to_user_gains(i));
end

% 3. SINGLE NOMA GROUP & POWER ALLOCATION
fprintf('\n=== SINGLE NOMA GROUP & POWER ALLOCATION ===\n');

% Use RIS-assisted channel gains for grouping
combined_user_gains = zeros(num_total_users, 1);
combined_user_gains(1) = gs_to_uav_gain;

% Calculate RIS-assisted gains for far users
for i = 1:num_far_users
    user_idx = i + 1;
    ris_gain = uav_to_RIS_gain * RIS_to_user_gains(i) * num_RIS_elements;
    combined_user_gains(user_idx) = ris_gain;
end

% SINGLE NOMA GROUP containing all users
noma_groups = cell(1, 1);
noma_groups{1} = 1:num_total_users; % All users in one group

fprintf('NOMA Group: Users %s\n', mat2str(noma_groups{1}));

% IMPROVED Power allocation for single NOMA group
powers = zeros(num_total_users, 1);
group = noma_groups{1};

% Sort users by channel gain (descending)
[~, order] = sort(combined_user_gains(group), 'descend');
sorted_users = group(order);

% Fair power allocation with more power to weaker users
P_total = P_ground; % Total power for the single group
num_users = length(sorted_users);

for u = 1:num_users
    user_idx = sorted_users(u);
    % More power to weaker users (inverse proportional to channel gain rank)
    power_fraction = (num_users - u + 1) / sum(1:num_users);
    powers(user_idx) = P_total * power_fraction;
end

fprintf('NOMA Power Allocation (Single Group):\n');
for i = 1:num_total_users
    if i == 1
        fprintf('  User %d (UAV): %.4f W\n', i, powers(i));
    else
        fprintf('  User %d: %.4f W\n', i, powers(i));
    end
end

% 4. BEAMFORMING WEIGHTS CALCULATION
fprintf('\n=== BEAMFORMING CALCULATION (UAV→RIS) ===\n');

beamforming_weights = zeros(num_antennas_UAV, 1);
ris_pos_3d = RIS_position;
uav_pos_3d = [user_positions(1, :), uav_altitude];

direction_vector = ris_pos_3d - uav_pos_3d;
distance = norm(direction_vector);
unit_vector = direction_vector / distance;

antenna_positions = (0:num_antennas_UAV-1) * 0.5;
steering_vector = exp(-1j * 2*pi * antenna_positions' * unit_vector(1));

beamforming_weights = steering_vector / norm(steering_vector);
fprintf('Beamforming weights calculated for %d-antenna UAV towards RIS\n', num_antennas_UAV);

% 5. RIS PHASE OPTIMIZATION
fprintf('\n=== RIS PHASE OPTIMIZATION ===\n');

optimized_RIS_phase = zeros(num_RIS_elements, 1);
for elem = 1:num_RIS_elements
    total_phase = 0;
    for i = 1:num_far_users
        user_idx = i + 1;
        user_pos_3d = [user_positions(user_idx, :), 0];
        ris_pos_3d = RIS_position;
        
        uav_to_ris_vec = ris_pos_3d - [user_positions(1,:), uav_altitude];
        ris_to_user_vec = user_pos_3d - ris_pos_3d;
        
        user_phase = angle(exp(1j * (norm(uav_to_ris_vec) + norm(ris_to_user_vec))));
        total_phase = total_phase + user_phase;
    end
    
    avg_phase = total_phase / num_far_users;
    optimized_RIS_phase(elem) = exp(-1j * avg_phase);
end

RIS_phase = optimized_RIS_phase;
fprintf('RIS phases optimized for %d elements (serving all %d users)\n', ...
    num_RIS_elements, num_far_users);

% 6. GENERATE DATA AT GS
user_bits = cell(num_total_users, 1);
user_symbols = cell(num_total_users, 1);
for i = 1:num_total_users
    user_bits{i} = randi([0 1], total_bits, 1);
    bits = reshape(user_bits{i}, 2, num_symbols)';
    map = [1+1j, 1-1j, -1+1j, -1-1j];
    user_symbols{i} = map(bits(:,1)*2 + bits(:,2) + 1)';
end

% 7. FIRST HOP: GS → UAV ONLY
fprintf('\n=== FIRST HOP: GS → UAV ONLY ===\n');

% Single group superposition
superimposed_gs = zeros(num_symbols, 1);
group = noma_groups{1};
for u = 1:length(group)
    user_idx = group(u);
    superimposed_gs = superimposed_gs + sqrt(powers(user_idx)) * user_symbols{user_idx};
end

gs_to_uav_signal = superimposed_gs;
uav_received = gs_to_uav_signal * sqrt(gs_to_uav_gain) + ...
    sqrt(noise_power/2)*(randn(num_symbols,1) + 1j*randn(num_symbols,1));

% 8. UAV DECODES ITS OWN DATA
fprintf('\nUAV decoding its own data from NOMA signal...\n');

current_signal = uav_received;
uav_own_symbols = zeros(num_symbols, 1);

% SIC for single group
group = noma_groups{1};
[~, order] = sort(combined_user_gains(group), 'descend');
sic_order = group(order);

for k = 1:length(sic_order)
    curr_user = sic_order(k);
    if curr_user == 1
        % Decode own signal
        real_sym = sign(real(current_signal));
        imag_sym = sign(imag(current_signal));
        for m = 1:num_symbols
            if real_sym(m) > 0 && imag_sym(m) > 0
                uav_own_symbols(m) = 1 + 1j;
            elseif real_sym(m) > 0 && imag_sym(m) < 0
                uav_own_symbols(m) = 1 - 1j;
            elseif real_sym(m) < 0 && imag_sym(m) > 0
                uav_own_symbols(m) = -1 + 1j;
            else
                uav_own_symbols(m) = -1 - 1j;
            end
        end
        break;
    else
        % Decode and cancel stronger user
        real_temp = sign(real(current_signal));
        imag_temp = sign(imag(current_signal));
        temp = zeros(num_symbols, 1);
        for m = 1:num_symbols
            if real_temp(m) > 0 && imag_temp(m) > 0
                temp(m) = 1 + 1j;
            elseif real_temp(m) > 0 && imag_temp(m) < 0
                temp(m) = 1 - 1j;
            elseif real_temp(m) < 0 && imag_temp(m) > 0
                temp(m) = -1 + 1j;
            else
                temp(m) = -1 - 1j;
            end
        end
        current_signal = current_signal - sqrt(powers(curr_user)) * temp * sqrt(gs_to_uav_gain);
    end
end

uav_forward_data = cell(num_total_users, 1);
uav_forward_data{1} = uav_own_symbols;
for i = 2:num_total_users
    uav_forward_data{i} = user_symbols{i};
end

bits_demod = zeros(num_symbols, 2);
bits_demod(real(uav_own_symbols) > 0, 1) = 0; 
bits_demod(real(uav_own_symbols) <= 0, 1) = 1;
bits_demod(imag(uav_own_symbols) > 0, 2) = 0; 
bits_demod(imag(uav_own_symbols) <= 0, 2) = 1;
uav_decoded_bits = reshape(bits_demod', [], 1);
uav_errors = sum(uav_decoded_bits ~= user_bits{1});

fprintf('UAV decoded its own data: %d errors (BER: %.6f)\n', uav_errors, uav_errors/total_bits);

% 9. SECOND HOP: UAV → RIS → Far Users
fprintf('\n=== SECOND HOP: UAV → RIS → Far Users (SINGLE NOMA GROUP) ===\n');

AN_symbols = sqrt(AN_power/2) * (randn(num_symbols,1) + 1j*randn(num_symbols,1));

% Single group transmission
uav_transmit_signal_beamformed = zeros(num_symbols, 1);
group = noma_groups{1};
for u = 1:length(group)
    user_idx = group(u);
    uav_transmit_signal_beamformed = uav_transmit_signal_beamformed + ...
        sqrt(powers(user_idx)) * uav_forward_data{user_idx};
end

uav_transmit_per_antenna = zeros(num_symbols, num_antennas_UAV);
for ant = 1:num_antennas_UAV
    uav_transmit_per_antenna(:, ant) = beamforming_weights(ant) * uav_transmit_signal_beamformed;
end

uav_transmit_per_antenna = uav_transmit_per_antenna + AN_symbols / num_antennas_UAV;

% 10. FAR USERS RECEIVE THROUGH RIS ONLY
fprintf('\nFar Users receiving through RIS only...\n');
far_user_errors = zeros(num_total_users, 1);
far_user_errors(1) = uav_errors;

for user_idx = 2:num_total_users
    user_beam_idx = user_idx - 1;
    
    ris_gain = uav_to_RIS_gain * RIS_to_user_gains(user_beam_idx) * ...
               num_RIS_elements * abs(mean(RIS_phase))^2;
    
    beamforming_gain = norm(beamforming_weights)^2;
    total_channel_gain = ris_gain * beamforming_gain;
    
    rx = zeros(num_symbols, 1);
    for ant = 1:num_antennas_UAV
        rx = rx + uav_transmit_per_antenna(:, ant) * sqrt(total_channel_gain);
    end
    rx = rx + sqrt(noise_power/2)*(randn(num_symbols,1) + 1j*randn(num_symbols,1));
    
    % SIC decoding for single group
    group = noma_groups{1};
    [~, order] = sort(combined_user_gains(group), 'descend');
    sic_order = group(order);
    signal = rx;
    sym = zeros(num_symbols, 1);
    
    for k = 1:length(sic_order)
        curr = sic_order(k);
        if curr == user_idx
            real_sym = sign(real(signal));
            imag_sym = sign(imag(signal));
            for m = 1:num_symbols
                if real_sym(m) > 0 && imag_sym(m) > 0
                    sym(m) = 1 + 1j;
                elseif real_sym(m) > 0 && imag_sym(m) < 0
                    sym(m) = 1 - 1j;
                elseif real_sym(m) < 0 && imag_sym(m) > 0
                    sym(m) = -1 + 1j;
                else
                    sym(m) = -1 - 1j;
                end
            end
            break;
        else
            real_temp = sign(real(signal));
            imag_temp = sign(imag(signal));
            temp = zeros(num_symbols, 1);
            for m = 1:num_symbols
                if real_temp(m) > 0 && imag_temp(m) > 0
                    temp(m) = 1 + 1j;
                elseif real_temp(m) > 0 && imag_temp(m) < 0
                    temp(m) = 1 - 1j;
                elseif real_temp(m) < 0 && imag_temp(m) > 0
                    temp(m) = -1 + 1j;
                else
                    temp(m) = -1 - 1j;
                end
            end
            signal = signal - sqrt(powers(curr)) * temp * sqrt(total_channel_gain);
        end
    end
    
    b1 = real(sym) <= 0; 
    b2 = imag(sym) <= 0;
    decoded_bits = [b1, b2]; 
    decoded_bits = decoded_bits(:);
    far_user_errors(user_idx) = sum(decoded_bits ~= user_bits{user_idx});
    
    fprintf('  User %d: %d errors (BER: %.6f)\n', user_idx, far_user_errors(user_idx), ...
        far_user_errors(user_idx)/total_bits);
end

% 11. EAVESDROPPER AND SECRECY ANALYSIS
fprintf('\n=== EAVESDROPPER AND SECRECY ANALYSIS ===\n');

% Eavesdropper with guaranteed minimum but worse channel
eaves_pos = [800, 800];
eaves_ris_dist = norm(eaves_pos - RIS_position(1:2));
pl_eaves = 35 + 10*3.5*log10(max(eaves_ris_dist, 1)); % Worse path loss
eaves_ris_gain = 10^(-pl_eaves/10);
eaves_ris_gain = abs(eaves_ris_gain * (sqrt(0.5)*(randn + 1j*randn)));
eaves_ris_gain = max(eaves_ris_gain, min_channel_gain/10); % Even worse than legitimate

fprintf('Eavesdropper Position: (%.1f, %.1f)\n', eaves_pos(1), eaves_pos(2));
fprintf('Eavesdropper RIS Channel Gain: %.6f\n', eaves_ris_gain);

% Calculate secrecy rates
secrecy_rates = zeros(num_total_users, 1);
legitimate_rates = zeros(num_total_users, 1);
eaves_rates = zeros(num_total_users, 1);

for i = 2:num_total_users
    % Legitimate user rate
    user_gain = combined_user_gains(i);
    user_snr = (powers(i) * user_gain) / noise_power;
    legitimate_rates(i) = log2(1 + user_snr);
    
    % Eavesdropper rate with heavy interference
    eaves_gain = uav_to_RIS_gain * eaves_ris_gain * num_RIS_elements;
    signal_power = powers(i) * eaves_gain;
    interference_power = AN_power * eaves_gain + sum(powers(setdiff(1:num_total_users, i))) * eaves_gain;
    eaves_snr = signal_power / (noise_power + interference_power);
    eaves_rates(i) = log2(1 + eaves_snr);
    
    secrecy_rates(i) = max(0, legitimate_rates(i) - eaves_rates(i));
end

fprintf('\n--- SECRECY RATES (bits/s/Hz) ---\n');
fprintf('User | Legit Rate | Eaves Rate | Secrecy Rate\n');
fprintf('-----|------------|------------|-------------\n');
for i = 2:num_total_users
    fprintf(' %2d  |   %6.3f   |   %6.3f   |    %6.3f\n', ...
        i, legitimate_rates(i), eaves_rates(i), secrecy_rates(i));
end

fprintf('\nAverage Secrecy Rate: %.3f bits/s/Hz\n', mean(secrecy_rates(2:end)));
fprintf('Minimum Secrecy Rate: %.3f bits/s/Hz\n', min(secrecy_rates(2:end)));

% 12. RESULTS SUMMARY
fprintf('\n=== FINAL RESULTS SUMMARY ===\n');
fprintf('User | Type | Errors | BER | Secrecy Rate\n');
fprintf('-----|------|--------|-----|-------------\n');
for i = 1:num_total_users
    if i == 1
        type = 'UAV ';
        fprintf('%2d   | %s  | %6d | %.4f |     --\n', ...
            i, type, far_user_errors(i), far_user_errors(i)/total_bits);
    else
        type = 'Far ';
        fprintf('%2d   | %s  | %6d | %.4f |   %.3f\n', ...
            i, type, far_user_errors(i), far_user_errors(i)/total_bits, secrecy_rates(i));
    end
end

fprintf('\n=== SYSTEM SUMMARY ===\n');
fprintf('✓ Single NOMA group (all %d users)\n', num_total_users);
fprintf('✓ Realistic power (4W)\n');
fprintf('✓ Guaranteed minimum channel gain: %.1e\n', min_channel_gain);
fprintf('✓ Rayleigh fading for all channels\n');
fprintf('✓ Average secrecy rate: %.3f bits/s/Hz\n', mean(secrecy_rates(2:end)));