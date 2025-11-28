%% 64-ELEMENT RIS + UAV + NOMA + AN (fixed)
clear; clc; close all; rng('shuffle');

fprintf('====================================================================\n');
fprintf('   64-ELEMENT RIS + UAV + NOMA + AN (fixed & robust)\n');
fprintf('====================================================================\n\n');

% ============ CHANGE THIS ONLY ============
N_users = 4;           % try 3,5,8,12,20
% =========================================

Ns      = 10000;
P_uav   = 10;
P_AN    = 1;
P_sig   = P_uav - P_AN;
N_ris   = 64;
noise   = 1e-12;

% Positions
uav_pos = [0, 0, 100];
ris_pos = [200, 0, 80];
theta   = linspace(0, 2*pi, N_users+1); theta(end) = [];
legit_pos = 150 * [cos(theta(1:N_users)); sin(theta(1:N_users)); zeros(1,N_users)].';
eve_pos = [500, 200, 0];

% Distances
d_uav_ris  = norm(uav_pos - ris_pos);
d_ris_user = sqrt(sum((legit_pos - ris_pos).^2, 2));
d_ris_eve  = norm(eve_pos - ris_pos);

% Channels (element-wise for RIS elements)
alpha_ur = 2.2; alpha_ru = 2.8;
h_ur = sqrt(d_uav_ris^(-alpha_ur)/2) .* (randn(N_ris,1) + 1j*randn(N_ris,1)); % N_ris x 1

h_ru = zeros(N_users, N_ris); % each row: 1 x N_ris
for u = 1:N_users
    h_ru(u,:) = sqrt(d_ris_user(u)^(-alpha_ru)/2) .* (randn(1, N_ris) + 1j*randn(1, N_ris));
end
h_re = sqrt(d_ris_eve^(-alpha_ru)/2) .* (randn(N_ris,1) + 1j*randn(N_ris,1)); % N_ris x 1

% RIS Phase design (coherent combination toward users average)
sum_h = sum(h_ru, 1);                     % 1 x N_ris
phase = exp(1j * angle(h_ur .* sum_h.')); % N_ris x 1
Theta = diag(phase);                      % N_ris x N_ris

% Effective channels (scalar complex) via RIS for each user and Eve
h_eff_legit = zeros(N_users,1);
for u = 1:N_users
    h_eff_legit(u) = h_ur.' * Theta * h_ru(u,:).';   % scalar (complex)
end
h_eff_eve = h_ur.' * Theta * h_re;                   % scalar complex

g_legit = abs(h_eff_legit).^2;
g_eve   = abs(h_eff_eve).^2;

% Print a quick sanity
fprintf('Eve effective gain: %.2e (%.2f dB)\n', g_eve, 10*log10(g_eve));
fprintf('Avg legit effective gain: %.2e (%.2f dB)\n\n', mean(g_legit), 10*log10(mean(g_legit)));

% Power allocation (weaker users get more)
[~, order] = sort(g_legit, 'ascend');
P = zeros(N_users,1);
for i = 1:N_users
    P(order(i)) = P_sig * (N_users + 1 - i) / sum(1:N_users);
end

% Generate bits and symbols (UAV + N_users)
total_users = N_users + 1;
bits_all = cell(total_users, 1);
syms_all = cell(total_users, 1);
qpsk = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);

for u = 1:total_users
    b = randi([0 1], 2*Ns, 1);
    bits_all{u} = b;
    % convert pairs to indices (left-msb)
    pairs = reshape(b, 2, Ns).';         % Ns x 2
    idx = 1 + pairs(:,1)*2 + pairs(:,2); % Ns x 1 (values 1..4)
    s = qpsk(idx).';                      % 1 x Ns -> transpose to get Ns x 1
    syms_all{u} = s(:);                   % ensure column vector (Ns x 1)
end

% Transmit (only ground users: syms_all{2}..syms_all{N_users+1})
tx = zeros(Ns, 1);
for u = 1:N_users
    tx = tx + sqrt(P(u)) * syms_all{u+1};  % Ns x 1
end
% Add AN (complex white) with specified power
tx = tx + sqrt(P_AN) * (randn(Ns,1) + 1j*randn(Ns,1));

% Receive: y_legit is N_users x Ns (each row is user's received vector)
% Use bsxfun for compatibility:
y_legit = bsxfun(@times, h_eff_legit, tx.');  % h_eff_legit: N_users x 1, tx.' 1 x Ns -> N_users x Ns
% add noise (complex)
y_legit = y_legit + sqrt(noise/2) * (randn(N_users, Ns) + 1j*randn(N_users, Ns));

% Eve receive (Ns x 1)
y_eve = (h_eff_eve .* tx) + sqrt(noise/2)*(randn(Ns,1) + 1j*randn(Ns,1));

% SIC + BER — fix: equalize before decision, ensure column vectors
ber = zeros(N_users,1);
for u = 1:N_users
    r_row = y_legit(u,:);         % 1 x Ns
    curr = r_row(:);              % Ns x 1 column working vector
    % Determine SIC decode order by received power at this user (P * |h|^2)
    recv_power = P .* (abs(h_eff_legit(u)).^2); % N_users x 1
    [~, sic] = sort(recv_power, 'descend');    % strongest-first
    
    for k = 1:N_users
        layer = sic(k);   % which user layer to decode now
        % Equalize by user's channel (complex scalar)
        rx_eq = curr ./ h_eff_legit(u);   % Ns x 1
        % Hard decision (normalized QPSK)
        dec = (sign(real(rx_eq)) + 1j*sign(imag(rx_eq))) / sqrt(2);
        dec = dec(:); % force column
        if layer == u
            % Demap bits from dec
            % Undo normalization for straightforward bit detection
            s = dec * sqrt(2);
            b1_hat = real(s) < 0;   % Ns x 1 logical
            b2_hat = imag(s) < 0;
            rec_bits = zeros(2*Ns,1);
            rec_bits(1:2:end) = b1_hat;
            rec_bits(2:2:end) = b2_hat;
            orig_bits = bits_all{u+1};
            len_min = min(length(rec_bits), length(orig_bits));
            errors = sum(rec_bits(1:len_min) ~= orig_bits(1:len_min));
            ber(u) = errors / len_min;
            break;
        else
            % Cancel this layer: subtract sqrt(P(layer)) * h * estimated symbol
            % Note: symbol estimate dec corresponds to current equalized estimate
            curr = curr - sqrt(P(layer)) * (h_eff_legit(u) .* dec);
        end
    end
end

% Secrecy rates (no AN vs with AN)
sec_no = zeros(N_users,1);
sec_yes = zeros(N_users,1);
for i = 1:N_users
    % interference at legit: sum of layers with smaller power (as per your original)
    I_legit = sum(P(P < P(i))) * g_legit(i);
    rate_l = log2(1 + P(i)*g_legit(i) / (I_legit + noise));
    % Eve without AN
    I_eve_no = (sum(P) - P(i)) * g_eve;
    rate_e_no = log2(1 + P(i)*g_eve / (I_eve_no + noise));
    % Eve with AN (AN adds P_AN * g_eve)
    rate_e_yes = log2(1 + P(i)*g_eve / (I_eve_no + P_AN*g_eve + noise));
    sec_no(i) = max(0, rate_l - rate_e_no);
    sec_yes(i) = max(0, rate_l - rate_e_yes);
end

% FINAL RESULT
fprintf('====================================================================\n');
fprintf('                        FINAL RESULT - %d USERS\n', N_users);
fprintf('====================================================================\n');
fprintf('Avg BER             → %.2e\n', mean(ber));
fprintf('Avg Secrecy (No AN) → %.6f bit/s/Hz\n', mean(sec_no));
fprintf('Avg Secrecy (RIS+AN)→ %.6f bit/s/Hz\n', mean(sec_yes));
if mean(sec_no) > 0
    fprintf('SECRECY GAIN        → +%.6f bit/s/Hz (+%.1f%%)\n', mean(sec_yes)-mean(sec_no), 100*(mean(sec_yes)/mean(sec_no)-1));
else
    fprintf('SECRECY GAIN        → +%.6f bit/s/Hz (baseline zero)\n', mean(sec_yes)-mean(sec_no));
end
fprintf('====================================================================\n\n');

% Print a short per-user summary (first 8 users)
fprintf('User | BER        | g_legit(dB) | dist_to_ris(m)\n');
for i = 1:min(8,N_users)
    fprintf('%4d | %.6e | %10.2f dB | %10.2f\n', i, ber(i), 10*log10(g_legit(i)), d_ris_user(i));
end
if N_users > 8, fprintf('... (total users = %d)\n', N_users); end
