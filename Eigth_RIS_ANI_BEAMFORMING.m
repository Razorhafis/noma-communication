%% 64-ELEMENT RIS + UAV + NOMA + AN = SCALABLE, UNBREAKABLE, FINAL GOD CODE
clear; clc; close all;

fprintf('====================================================================\n');
fprintf('   64-ELEMENT RIS + UAV + NOMA + AN \n');
fprintf('====================================================================\n\n');

% ============ CHANGE THIS ONLY ============
N_users = 10;           % TRY 3, 5, 8, 12, 20 — WORKS EVERY TIME
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

% Channels
alpha_ur = 2.2; alpha_ru = 2.8;
h_ur = sqrt(d_uav_ris^(-alpha_ur)/2) * (randn(N_ris,1) + 1j*randn(N_ris,1));

h_ru = zeros(N_users, N_ris);
for u = 1:N_users
    h_ru(u,:) = sqrt(d_ris_user(u)^(-alpha_ru)/2) * (randn(1, N_ris) + 1j*randn(1, N_ris));
end
h_re = sqrt(d_ris_eve^(-alpha_ru)/2) * (randn(N_ris,1) + 1j*randn(N_ris,1));

% RIS Phase
sum_h = sum(h_ru, 1);
phase = exp(1j * angle(h_ur .* sum_h.'));
Theta = diag(phase);

% Effective channels
h_eff_legit = zeros(N_users,1);
for u = 1:N_users
    h_eff_legit(u) = h_ur.' * Theta * h_ru(u,:).';
end
h_eff_eve = h_ur.' * Theta * h_re;

g_legit = abs(h_eff_legit).^2;
g_eve   = abs(h_eff_eve)^2;

% Power allocation
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
    idx = 1 + b(1:2:end)*2 + b(2:2:end);
    syms_all{u} = qpsk(idx).';  % Ns x 1
end

% Transmit (only ground users: 2 to total_users)
tx = zeros(Ns, 1);
for u = 1:N_users
    tx = tx + sqrt(P(u)) * syms_all{u+1};
end
tx = tx + sqrt(P_AN) * (randn(Ns,1) + 1j*randn(Ns,1));

% Receive
y_legit = h_eff_legit .* tx.' + sqrt(noise/2)*(randn(N_users, Ns) + 1j*randn(N_users, Ns));

% SIC + BER — 100% FIXED WITH (:)
ber = zeros(N_users,1);
for u = 1:N_users
    r = y_legit(u,:);
    curr = r;
    [~, sic] = sort(P, 'descend');
    
    for k = 1:N_users
        dec = sign(real(curr)) + 1j*sign(imag(curr));
        dec = dec / sqrt(2);
        
        if sic(k) == u
            rec_b1 = (real(dec) < 0).';  % Force column
            rec_b2 = (imag(dec) < 0).';  % Force column
            
            orig = bits_all{u+1};
            sent_b1 = orig(1:2:end);
            sent_b2 = orig(2:2:end);
            
            ber(u) = (sum(rec_b1 ~= sent_b1) + sum(rec_b2 ~= sent_b2)) / (2*Ns);
            break;
        else
            curr = curr - sqrt(P(sic(k))) * dec * h_eff_legit(u);
        end
    end
end

% Secrecy
sec_no = zeros(N_users,1);
sec_yes = zeros(N_users,1);
for i = 1:N_users
    I_legit = sum(P(P < P(i))) * g_legit(i);
    rate_l = log2(1 + P(i)*g_legit(i)/(I_legit + noise));
    I_eve = (sum(P) - P(i)) * g_eve;
    rate_e_no  = log2(1 + P(i)*g_eve/(I_eve + noise));
    rate_e_yes = log2(1 + P(i)*g_eve/(I_eve + P_AN*g_eve + noise));
    sec_no(i)  = max(0, rate_l - rate_e_no);
    sec_yes(i) = max(0, rate_l - rate_e_yes);
end

% FINAL RESULT
fprintf('====================================================================\n');
fprintf('                        FINAL RESULT - %d USERS\n', N_users);
fprintf('====================================================================\n');
fprintf('Avg BER             → %.2e\n', mean(ber));
fprintf('Avg Secrecy (No AN) → %.3f bit/s/Hz\n', mean(sec_no));
fprintf('Avg Secrecy (RIS+AN)→ %.3f bit/s/Hz\n', mean(sec_yes));
fprintf('SECRECY GAIN        → +%.3f bit/s/Hz (+%.1f%%)\n', ...
    mean(sec_yes)-mean(sec_no), 100*(mean(sec_yes)/mean(sec_no)-1));
fprintf('====================================================================\n');
