clear all;
close all;

% Monte Carlo Parameters
num_monte_carlo = 50;
mc_results.sinr = zeros(num_monte_carlo, 10);
mc_results.rate = zeros(num_monte_carlo, 10);
mc_results.secrecy_rate = zeros(num_monte_carlo, 10);
mc_results.powers = zeros(num_monte_carlo, 10);
mc_results.channel_gains = zeros(num_monte_carlo, 10);

% System Parameters
num_users = 10;
P_total = 2;
AN_power_fraction = 0.03;
AN_power = AN_power_fraction * P_total;
noise_power = 1e-6;

% NOMA Parameters
num_noma_groups = 5;
users_per_group = 2;

% Monte Carlo Simulation Loop
for mc_iter = 1:num_monte_carlo
    fprintf('NOMA Iteration %d/%d\n', mc_iter, num_monte_carlo);
    
    % 1. Generate random channel gains
    channel_gains = 0.1 + 0.9 * rand(num_users, 1);
    mc_results.channel_gains(mc_iter, :) = channel_gains';
    
    % 2. NOMA User Grouping
    [~, user_ordered_indices] = sort(channel_gains, 'descend');
    noma_groups = cell(num_noma_groups, 1);
    for g = 1:num_noma_groups
        start_idx = (g-1)*users_per_group + 1;
        end_idx = min(g*users_per_group, num_users);
        noma_groups{g} = user_ordered_indices(start_idx:end_idx);
    end
    
    % 3. NOMA Power Allocation
    powers = zeros(num_users, 1);
    for g = 1:num_noma_groups
        group_users = noma_groups{g};
        if length(group_users) == 2
            strong_user = group_users(1);
            weak_user = group_users(2);
            powers(weak_user) = 0.7 * (P_total - AN_power) / num_noma_groups;
            powers(strong_user) = 0.3 * (P_total - AN_power) / num_noma_groups;
        else
            powers(group_users) = (P_total - AN_power) / num_noma_groups;
        end
    end
    mc_results.powers(mc_iter, :) = powers';
    
    % 4. Calculate rates for each user
    for user_idx = 1:num_users
        % Find user's group
        user_group = 0;
        for g = 1:num_noma_groups
            if ismember(user_idx, noma_groups{g})
                user_group = g;
                break;
            end
        end
        
        if user_group == 0
            continue;
        end
        
        % LEGITIMATE USER RATE
        signal_power = powers(user_idx) * channel_gains(user_idx)^2;
        sinr_main = signal_power / noise_power;
        rate_main = log2(1 + sinr_main);
        
        % EAVESDROPPER RATE
        other_users = setdiff(1:num_users, user_idx);
        [~, eaves_idx] = max(channel_gains(other_users));
        eavesdropper_gain = channel_gains(other_users(eaves_idx));
        
        % Eavesdropper experiences NOMA interference + AN
        noma_interference = 0;
        for g = 1:num_noma_groups
            if g ~= user_group
                group_power = sum(powers(noma_groups{g}));
                noma_interference = noma_interference + group_power * eavesdropper_gain^2;
            end
        end
        
        an_interference = AN_power * eavesdropper_gain^2;
        total_interference = noma_interference + an_interference + noise_power;
        
        sinr_eaves = (powers(user_idx) * eavesdropper_gain^2) / total_interference;
        rate_eaves = log2(1 + sinr_eaves);
        
        % SECRECY RATE
        secrecy_rate = max(0, rate_main - rate_eaves);
        
        % Store results
        mc_results.sinr(mc_iter, user_idx) = 10 * log10(sinr_main);
        mc_results.rate(mc_iter, user_idx) = rate_main;
        mc_results.secrecy_rate(mc_iter, user_idx) = secrecy_rate;
    end
end

% =========================================================================
% DISPLAY RESULTS
% =========================================================================

fprintf('\n\n');
fprintf('===================================================\n');
fprintf('      NOMA SECRECY SIMULATION RESULTS\n');
fprintf('===================================================\n\n');

% Display detailed results for each user
fprintf('=== PER-USER PERFORMANCE ===\n');
fprintf('User | Avg SINR | Avg Rate | Sec Rate |  Power  | Chan Gain\n');
fprintf('-----|----------|----------|----------|---------|----------\n');

for user_idx = 1:num_users
    avg_sinr = mean(mc_results.sinr(:, user_idx));
    avg_rate = mean(mc_results.rate(:, user_idx));
    avg_secrecy = mean(mc_results.secrecy_rate(:, user_idx));
    avg_power = mean(mc_results.powers(:, user_idx));
    avg_gain = mean(mc_results.channel_gains(:, user_idx));
    
    fprintf(' %2d  |  %6.2f  |  %6.3f  |  %6.3f  | %7.4f |  %6.4f\n', ...
        user_idx, avg_sinr, avg_rate, avg_secrecy, avg_power, avg_gain);
end

% Strong vs Weak user analysis
fprintf('\n=== NOMA GROUP ANALYSIS ===\n');
strong_users = [1, 3, 5, 7, 9];
weak_users = [2, 4, 6, 8, 10];

strong_sinr = mean(mean(mc_results.sinr(:, strong_users)));
strong_secrecy = mean(mean(mc_results.secrecy_rate(:, strong_users)));
strong_power = mean(mean(mc_results.powers(:, strong_users)));

weak_sinr = mean(mean(mc_results.sinr(:, weak_users)));
weak_secrecy = mean(mean(mc_results.secrecy_rate(:, weak_users)));
weak_power = mean(mean(mc_results.powers(:, weak_users)));

fprintf('Strong Users (Better channels, Less power):\n');
fprintf('  Avg SINR: %.2f dB, Avg Secrecy: %.3f bps/Hz, Avg Power: %.4f W\n', ...
    strong_sinr, strong_secrecy, strong_power);

fprintf('Weak Users (Worse channels, More power):\n');
fprintf('  Avg SINR: %.2f dB, Avg Secrecy: %.3f bps/Hz, Avg Power: %.4f W\n', ...
    weak_sinr, weak_secrecy, weak_power);

% System statistics
fprintf('\n=== SYSTEM STATISTICS ===\n');
fprintf('Total Average Rate:        %.3f bps/Hz\n', mean(sum(mc_results.rate, 2)));
fprintf('Total Average Secrecy Rate: %.3f bps/Hz\n', mean(sum(mc_results.secrecy_rate, 2)));
fprintf('Secrecy Efficiency:        %.1f%%\n', ...
    (mean(sum(mc_results.secrecy_rate, 2)) / mean(sum(mc_results.rate, 2))) * 100);

positive_secrecy = sum(mc_results.secrecy_rate(:) > 0) / numel(mc_results.secrecy_rate) * 100;
fprintf('Cases with Positive Secrecy: %.1f%%\n', positive_secrecy);

% =========================================================================
% VISUALIZATION (FIXED - NO ALPHA PROPERTY)
% =========================================================================

figure('Position', [100, 100, 1200, 800]);

% 1. SINR comparison
subplot(2,3,1);
avg_sinr_per_user = mean(mc_results.sinr, 1);
bar(avg_sinr_per_user);
title('Average SINR per User');
xlabel('User Index');
ylabel('SINR (dB)');
grid on;

% 2. Secrecy rate comparison
subplot(2,3,2);
avg_secrecy_per_user = mean(mc_results.secrecy_rate, 1);
bar(avg_secrecy_per_user);
title('Average Secrecy Rate per User');
xlabel('User Index');
ylabel('Secrecy Rate (bps/Hz)');
grid on;

% 3. Power allocation
subplot(2,3,3);
avg_power_per_user = mean(mc_results.powers, 1);
bar(avg_power_per_user);
title('Power Allocation per User');
xlabel('User Index');
ylabel('Power (W)');
grid on;

% 4. Channel gains vs Power allocation
subplot(2,3,4);
avg_gains = mean(mc_results.channel_gains, 1);
scatter(avg_gains, avg_power_per_user, 100, 'filled');
for i = 1:num_users
    text(avg_gains(i), avg_power_per_user(i), sprintf('U%d', i), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
xlabel('Channel Gain');
ylabel('Power (W)');
title('NOMA Principle: More Power to Weaker Channels');
grid on;

% 5. SINR vs Secrecy Rate (FIXED LINE - no alpha)
subplot(2,3,5);
scatter(mc_results.sinr(:), mc_results.secrecy_rate(:), 30, 'filled');
xlabel('SINR (dB)');
ylabel('Secrecy Rate (bps/Hz)');
title('SINR vs Secrecy Rate');
grid on;

% 6. Rate comparison
subplot(2,3,6);
avg_rates = mean(mc_results.rate, 1);
bar_data = [avg_rates; avg_secrecy_per_user]';
bar(bar_data);
legend('Achievable Rate', 'Secrecy Rate', 'Location', 'best');
title('Rate Comparison per User');
xlabel('User Index');
ylabel('Rate (bps/Hz)');
grid on;

fprintf('\n=== SIMULATION SUMMARY ===\n');
fprintf('This demonstrates a secure NOMA system where:\n');
fprintf('- Users are grouped by channel quality\n');
fprintf('- Weak users get MORE power (fairness)\n');
fprintf('- Artificial Noise + NOMA interference provides security\n');
fprintf('- Secrecy rates show protected communication\n');

fprintf('\nSimulation completed successfully!\n');