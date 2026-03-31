% MAE 283A Final Project
% Charlie Frank
% Fall 2025
% Load the .mat file
clear
close all
clc
load('step_Pinv_Pbess.mat');

%% Parametric Identification Setup 
    t = t(:); % TIME
    u = Pin_INV; % INPUT
    y = Pout_BESS; % OUTPUT 
    N = length(u); % SIZE
    Ts = median(diff(t)); % T_sample
    Fs = 1/Ts; % F_sample

% Fig 1 values for plot
    [Ru_hat, lags] = xcorr(u, 'biased');
    tau = lags * Ts;

% Fig 2 values for plot
    u = u - mean(u);
    U = fft(u);
    Phi_u = (Ts/N) * abs(U(1:floor(N/2)+1)).^2;
    w = 2*pi*(0:floor(N/2))'/N*Fs;

% Fig 3
    y = y - mean(y);
    Y = fft(y);
    Phi_y  = (Ts/N) * abs(Y(1:floor(N/2)+1)).^2;
    Phi_yu = (Ts/N) * (Y(1:floor(N/2)+1) .* conj(U(1:floor(N/2)+1))); 
    Ghat = Phi_yu ./ Phi_u;
    
    w_min = 0.025;
    w_max = 2.38745;
    idx = (w >= w_min) & (w <= w_max);
% Fig 4
    n = 435;
    % Final FIR estimation with chosen order
    Phi = zeros(N-n, n+1);
    col = u(n+1:N);
    row = u(n+1:-1:1);
    Phi = toeplitz(col, row);

    y_vec = y(n+1:N);
    g_hat = Phi \ y_vec;
    residuals = y_vec - Phi*g_hat;
    sigma2 = var(residuals);
    SE_full = sqrt(sigma2 * diag(inv(Phi'*Phi)));
    alpha = 2.576 * SE_full; % 99% confidence bounds

    k_vec = 0:n;
% Fig 5
    eps_hat = residuals; 
    u_trunc = u(n+1:end);
    [Reu_hat, lags] = xcorr(eps_hat, u_trunc, 2*n, 'biased'); 
    tau2 = lags * Ts;
    Reu_hat = Reu_hat / (std(eps_hat) * std(u));
    conf99 = 2.576 / sqrt(N); 
    idx_pos = lags >= 0;
    tau2_pos = tau2(idx_pos);
    Reu_hat_pos = Reu_hat(idx_pos);
    % Delay
    lags_pos = lags(idx_pos);
    Re_pos = Reu_hat(idx_pos);

    hit = find(abs(Re_pos) > conf99, 1, 'first');
    if ~isempty(hit)
        k_delay = lags_pos(hit);        % in samples
        tau_delay = k_delay * Ts;       % in seconds
    else
        k_delay = 0;
        tau_delay = 0;
    end
    fprintf('Estimated delay: %d samples (%.3f s)\n', k_delay, tau_delay);
% Fig 6
    g_hat = g_hat(:);
    D = g_hat(1);
    s = min(ceil(n/6), 200);   % number of block rows
    p = min(ceil(n/6), 200);   % number of block columns

    H0 = zeros(s, p);
    H1 = zeros(s, p);
    for i = 1:s
        for j = 1:p
            H0(i,j) = g_hat(i + j - 1);
            H1(i,j) = g_hat(i + j);
        end
    end
    % SVD of H0
    [U, S, V] = svd(H0, 'econ');
    sigma = diag(S);

% Fig 7
thr = 0.2248*max(sigma) ;              % relative threshold
r = max(1, sum(sigma > thr));         % ensure at least order-1
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
Sr_sqrt= sqrt(Sr);
Sr_inv_sqrt = diag(1./diag(Sr_sqrt));
A = Sr_inv_sqrt * (Ur' * H1 * Vr) * Sr_inv_sqrt;
B = Sr_sqrt * Vr(1,:).';
C = (Ur(1,:) * Sr_sqrt);
sysd = ss(A, B, C, D, Ts);
Kmax = 200;   % number of samples to display
g_est = zeros(Kmax,1);
g_est(1) = D;   % k=0 term
for k = 2:Kmax
    g_est(k) = C * (A^(k-1)) * B;
end

% Fig 8
    y_ss = lsim(sysd, u, t);   % state-space simulated output
    eps_ss = y - y_ss;         % simulation error
    [Reu_ss_hat, lags] = xcorr(eps_ss, u, 2*r, 'unbiased');
    tau3 = lags * Ts;
    idx_pos2 = lags >= 0;
    tau3_pos = tau3(idx_pos2);
    Reu_ss_hat = Reu_ss_hat(idx_pos2);
    % Normalize cross-correlation
    Reu_ss_hat = Reu_ss_hat / (std(eps_ss) * std(u));

    % 99% confidence bounds
    conf99 = 2.576 / sqrt(N);

% Fig 9 >> Box Jenkins model various values tested 33330 found to be best
    data = iddata(y, u, Ts);
    nb = 3; nc = 3; nd = 3; nf = 3; nk = 0;
    modelBJ = bj(data, [nb nc nd nf nk]);
    eps_id = pe(data, modelBJ);
    epsilon = eps_id.OutputData;
    N = length(epsilon);
    eps_zm = epsilon - mean(epsilon);
    u_zm   = u - mean(u);
    maxLag = 40;
    [Reu, lags] = xcorr(eps_zm, u_zm, maxLag, 'unbiased');
    rho = Reu ./ (std(eps_zm)*std(u_zm) + eps);
    CI = 2.58/sqrt(N);

% Fig 10
    y_sim = sim(modelBJ, u);
    y_pred_id = predict(modelBJ, data, 1);
    y_pred = y_pred_id.OutputData;

%% FIGURE 1 and Explanation:
figure(1);
plot(tau, Ru_hat, 'b-', 'LineWidth', 1.2);
grid on;
xlabel('\tau (s)');
ylabel('hat{R}^N_u(\tau)');
title('Figure 1: Estimate hat{R}^N_u(\tau) of Input Autocorrelation Function R_u(\tau)');

%% FIGURE 2 and Explanation:
figure(2); 
plot(w, (Phi_u), 'LineWidth', 1.2);
xlabel('\omega  [rad/s]'); 
ylabel('hat{\Phi}^N_u(\omega)');
title('Estimated Input Spectrum');
grid on;

% 0.025 to 2.38745 rad/s reliable range visual inspection
 %% FIGURE 3 and Explanation:
figure(3);
subplot(2,1,1)
 %hold on
semilogx(w(idx), 20*log10(abs(Ghat(idx))), 'LineWidth', 1.2);
       % semilogx(w(idx), 20*log10(abs(Ghat_smooth(idx))), 'r','LineWidth', 1.2);
ylabel('Magnitude (dB)');
title('SPA Estimate of G(\omega)');
grid on;

subplot(2,1,2)
 %hold on
semilogx(w(idx), unwrap(angle(Ghat(idx))) * 180/pi, 'LineWidth', 1.2);
      %  semilogx(w(idx), unwrap(angle(Ghat_smooth(idx)))*180/pi,'r', 'LineWidth', 1.2);
xlabel('\omega  [rad/s]');
ylabel('Phase (deg)');
grid on;

%% FIGURE 4 and Explanation
figure(4);
hold on;
plot(k_vec, g_hat, 'b--', 'LineWidth', 1.5, 'DisplayName', 'g_N(k)');
plot(k_vec, g_hat + alpha, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Upper 99% bound');
plot(k_vec, g_hat - alpha, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Lower 99% bound');
xlabel('k'); 
ylabel('g_N(k)');
title(['Figure 4: Estimated FIR Impulse Response (n = ' num2str(n) ') with 99% Confidence Bounds']);
grid on;
legend('Location','best');

%% Figure 5 and Explanation
figure(5); 
stem(tau2_pos, Reu_hat_pos, 'b--', 'filled'); hold on;
yline(conf99, 'r--', 'LineWidth', 1.2, 'DisplayName','99% bound');
yline(-conf99, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('\tau (s)');
ylabel('hat{R}^N_{\epsilon u}(\tau)');
title(['Figure 5: Cross-correlation R^N_{\epsilon u}(\tau) with 99% Confidence Bounds, n = ' num2str(n)]);
legend('Cross-corr','99% bounds','Location','best');
% Generally within bounds, then exits range arou

%% FIGURE 6 and Explanation:
figure(6);
semilogy(1:length(sigma), sigma, 'bo-','LineWidth',1.5,'MarkerSize',5);
grid on;
xlabel('Order index');
ylabel('Hankel singular value');
title('Figure 6: ERA-based Hankel singular values');

%% FIGURE 7 and Explanation:
figure(7);
stem(0:Kmax-1, g_est, 'b','filled');
grid on;
xlabel('k');
ylabel('g(k)');
title(sprintf('Figure 7: Impulse Response of Estimated State-Space Model (order = %d)', size(A,1)));

%% FIGURE 8 and Explanation:
figure(8);
stem(tau3_pos, Reu_ss_hat, 'b--','filled'); hold on;
yline(conf99, 'r--','LineWidth',1.2,'DisplayName','99% bound');
yline(-conf99, 'r--','LineWidth',1.2);
grid on;
xlabel('\tau (s)');
ylabel('\hat{R}^N_{\epsilon_{SS}u}(\tau)');
title(['Figure 8: Cross-correlation R^N_{\epsilon_{SS}u}(\tau) with 99% Confidence Bounds, n = ' num2str(r)]);
legend('Cross-corr','99% bounds','Location','best');

%% FIGURE 9 and Explanation:
figure(9); 
stem(lags, rho, 'filled'); hold on;
yline(+CI, 'r--', '99% CI'); yline(-CI, 'r--', '99% CI');
xline(+2*max(nb,nf), 'k:'); xline(-2*max(nb,nf), 'k:');
xlabel('Lag τ'); ylabel('ρ_{εu}(τ)');
title('Figure 9: Residual–input correlation (BJ [3 3 3 3 0])');
grid on;

%% FIGURE 10 and Explanation:
figure(10); 
plot(t, y, 'k-', 'LineWidth', 1.2); hold on;
plot(t, y_sim, 'b--', 'LineWidth', 1.2);
plot(t, y_pred, 'r-', 'LineWidth', 1.2);
legend('Measured y(t)', 'Simulated y_{sim}(t)', 'Predicted y(t|t-1)', 'Location', 'best');
xlabel('Time [s]'); ylabel('Power');
title('Figure 10: Measured vs Simulation vs Prediction (BJ [3 3 3 3 0])');
grid on;