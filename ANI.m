% TWO-HOP NOMA + DIRECTIONAL AN (NULL STEERING) - FIXED
clear; close all; clc; rng(42);  % FIXED SEED FOR CONSISTENT RESULTS

fprintf('===== DIRECTIONAL AN NULL STEERING (fixed) =====\n\n');

%% Parameters
N_users = 10;        % number of legitimate ground users
Ns = 10000;          % number of symbols
P_uav = 10;
AN_frac = 0.05;
P_AN = AN_frac * P_uav;
P_sig = P_uav - P_AN;
noise = 1e-10;

%% Geometry (UAV at origin projection) - FIXED POSITIONS
uav_pos = [0,0,100];

% FIXED user positions for consistent results
fixed_angles = [0, 45, 90, 135, 180, 225, 270, 315, 30, 60] * pi/180;
fixed_radii = [120, 110, 130, 115, 125, 105, 135, 100, 140, 95];

legit_pos = zeros(N_users, 2);
for u = 1:N_users
    angle = fixed_angles(u);
    radius = fixed_radii(u);
    legit_pos(u,:) = [radius * cos(angle), radius * sin(angle)];
end

eve_pos = [350, 100, 0];

% Distances (2D distances from UAV projection)
d_legit = sqrt(sum((legit_pos - [0,0]).^2, 2));   % N_users x 1
d_eve = norm(eve_pos(1:2) - [0,0]);               % scalar

% Path loss
alpha_legit = 3.2;
alpha_eve = 4.0;
PL_legit = d_legit.^(-alpha_legit);
PL_eve = d_eve^(-alpha_eve);

% Small-scale Rayleigh fading (complex channel coefficients)
h_legit = sqrt(PL_legit/2) .* (randn(N_users,1) + 1j*randn(N_users,1)); % N_users x 1
g_legit = abs(h_legit).^2;
h_eve = sqrt(PL_eve/2) * (randn + 1j*randn); % scalar
g_eve = abs(h_eve)^2;

fprintf('Eve gain: %.2e (%.1f dB)\n', g_eve, 10*log10(g_eve));
fprintf('Avg legit gain: %.2e (%.1f dB)\n\n', mean(g_legit), 10*log10(mean(g_legit)));

%% Power allocation (weaker users get more power)
[~, order] = sort(g_legit, 'ascend');
P = zeros(N_users,1);
for i = 1:N_users
    P(order(i)) = P_sig * (N_users + 1 - i) / sum(1:N_users);
end

%% Generate QPSK symbols (no comm toolbox)
qpsk_const = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2); % normalized
bits = cell(N_users+1,1);
syms = cell(N_users+1,1);

for u = 1:(N_users+1)
    b = randi([0 1], Ns*2, 1);
    bits{u} = b;
    % reshape to (Ns x 2) and convert to indices 1..4 (left-msb)
    idx = 1 + bi2de(reshape(b, 2, Ns).', 'left-msb'); % 1 x Ns
    s = qpsk_const(idx).'; % column vector Ns x 1
    syms{u} = s(:);        % force column
end

%% DIRECTIONAL AN: Project INTO NULL SPACE OF LEGITIMATE USERS
% Unit vectors from UAV projection to users (N_users x 2)
directions = legit_pos ./ d_legit; % elementwise (N_users x 2)
H = directions.';                  % 2 x N_users

% SVD-based null direction (gives 2x2 U, take last column)
[U,S,V] = svd(H,'econ');            % U: 2x2, V: N_users x N_users (if N_users>2)
null_space = U(:, end);             % 2 x 1
P_null = null_space * null_space.'; % 2 x 2 projection
if trace(P_null) == 0
    P_null = eye(2); % fallback
end
P_null = P_null / trace(P_null);    % normalize

% Generate AN sequence and project
z = (randn(Ns,1) + 1j*randn(Ns,1)) * sqrt(P_AN); % complex Ns x 1
an_2d = [real(z), imag(z)];                      % Ns x 2
an_projected = (P_null * an_2d.').';             % Ns x 2
an_signal = an_projected(:,1) + 1j*an_projected(:,2); % Ns x 1 complex

%% Build transmit signal (superposition of user symbols) - column vector
tx = zeros(Ns,1);
for i = 1:N_users
    tx = tx + sqrt(P(i)) * syms{i+1};  % syms{2..N_users+1} are user symbols
end
tx = tx + an_signal; % add AN

%% Receive at legitimate users and Eve
% We'll store received signals as Nx: Ns x N_users (columns)
rx_legit = zeros(Ns, N_users);  % columns: each legitimate user's received vector
for i = 1:N_users
    rx_legit(:,i) = tx .* h_legit(i) + sqrt(noise/2) * (randn(Ns,1) + 1j*randn(Ns,1));
end
rx_eve = tx * h_eve + sqrt(noise/2)*(randn(Ns,1) + 1j*randn(Ns,1)); % Ns x 1

%% SIC Decoding at legitimate users (equalize -> decode -> cancel)
ber_legit = zeros(N_users,1);
for target = 1:N_users
    rx = rx_legit(:,target);      % Ns x 1
    cur = rx;                     % working copy
    % Determine SIC decode order by received power at this user (descending)
    recv_power_layers = P .* (abs(h_legit(target))^2);  % N_users x 1
    [~, sic_order] = sort(recv_power_layers, 'descend'); % strongest-first
    
    for k = 1:N_users
        j = sic_order(k);  % layer index (1..N_users)
        
        % Equalize by user's channel before hard decision
        rx_eq = cur ./ h_legit(target);   % Ns x 1
        sym_hat = (sign(real(rx_eq)) + 1j * sign(imag(rx_eq))) / sqrt(2); % Ns x 1
        sym_hat = sym_hat(:);
        
        if j == target
            % Demap to bits
            s = sym_hat * sqrt(2); % undo normalization for bit test
            b1 = real(s) < 0;
            b2 = imag(s) < 0;
            rec_bits = zeros(Ns*2,1);
            rec_bits(1:2:end) = b1;
            rec_bits(2:2:end) = b2;
            % Compare safely
            tx_bits = bits{target+1};
            len_min = min(length(rec_bits), length(tx_bits));
            ber_legit(target) = mean(rec_bits(1:len_min) ~= tx_bits(1:len_min));
            break;
        else
            % Cancel this layer: subtract sqrt(P(j)) * h * s_hat
            cur = cur - sqrt(P(j)) * (h_legit(target) .* sym_hat);
        end
    end
end

%% Secrecy Rate Calculation (with and without AN)
sec_no_an = zeros(N_users,1);
sec_with_an = zeros(N_users,1);

for i = 1:N_users
    % interference at legit: sum of lower-power layers (as in your original code)
    interference_legit = 0;
    for j = 1:N_users
        if j ~= i && P(j) < P(i)
            interference_legit = interference_legit + P(j) * g_legit(i);
        end
    end
    sig_leg = P(i) * g_legit(i);
    snr_leg_no = sig_leg / (interference_legit + noise);
    rate_leg_no = log2(1 + snr_leg_no);
    
    % eavesdropper without AN
    interference_eve_no = (sum(P) - P(i)) * g_eve;
    snr_eve_no = (P(i) * g_eve) / (interference_eve_no + noise);
    rate_eve_no = log2(1 + snr_eve_no);
    sec_no_an(i) = max(0, rate_leg_no - rate_eve_no);
    
    % With AN: assume AN is nulled at legit users (ideal null steering), so legit unaffected.
    snr_leg_an = sig_leg / (interference_legit + 1e-12 + noise); % tiny floor to avoid /0
    rate_leg_an = log2(1 + snr_leg_an);
    % Eve sees AN power scaled by its channel power: P_AN * g_eve
    interference_eve_an = interference_eve_no + P_AN * g_eve;
    snr_eve_an = (P(i) * g_eve) / (interference_eve_an + noise);
    rate_eve_an = log2(1 + snr_eve_an);
    sec_with_an(i) = max(0, rate_leg_an - rate_eve_an);
end

%% Final results
fprintf('====================================================\n');
fprintf('               FINAL RESULT (fixed run)\n');
fprintf('====================================================\n');
fprintf('Avg Secrecy (No AN)   : %.6f bits/s/Hz\n', mean(sec_no_an));
fprintf('Avg Secrecy (With AN) : %.6f bits/s/Hz\n', mean(sec_with_an));
improvement = mean(sec_with_an) - mean(sec_no_an);
if mean(sec_no_an) > 0
    pct = 100 * (mean(sec_with_an)/mean(sec_no_an) - 1);
else
    pct = Inf;
end
fprintf('IMPROVEMENT: +%.6f bits/s/Hz (+%.1f%%)\n', improvement, pct);
fprintf('Avg User BER: %.6e\n', mean(ber_legit));
fprintf('====================================================\n');

% Optional: print per-user summary (first few users to keep terminal readable)
fprintf('\nUser | BER        | g(dB)    | dist(m) \n');
for i = 1:min(8, N_users)
    fprintf('%4d | %.6e | %7.2f dB | %7.2f\n', i, ber_legit(i), 10*log10(g_legit(i)), d_legit(i));
end
fprintf('... (total users = %d)\n\n', N_users);
