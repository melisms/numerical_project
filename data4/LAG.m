%% =========================================================
%  LAG_data4.m
%  En iyi zaman offsetini bul, sonra PLS uygula — data4
%  (Hic Toolbox gerektirmez)
% =========================================================
clearvars
clc

%% =========================
% 1. OZELLIK VERiSiNi YUKLE
% =========================
load('3rd_fit_data4.mat', 'T');
X = table2array(T);

%% =========================
% 2. AKIS VERiSiNi YUKLE
% =========================
pp = load('pp.mat');
sweep_time_s = pp.R_time(:) * 60;   % dakika -> saniye

flow_path = 'C:\Users\asli\Desktop\numerical_project-main\IWIS2025flow\Auto_saved_file4.txt';
opts = detectImportOptions(flow_path, 'Delimiter', '\t');
opts.VariableNamingRule = 'modify';
opts.Encoding = 'ISO-8859-1';
Traw = readtable(flow_path, opts);

flow_t  = Traw{:,2};
s1_tgt  = Traw{:,3};
s2_tgt  = Traw{:,5};

ok = ~isnan(flow_t) & ~isnan(s1_tgt) & ~isnan(s2_tgt);
flow_t = flow_t(ok);
s1_tgt = s1_tgt(ok);
s2_tgt = s2_tgt(ok);

total_flow = s1_tgt + s2_tgt;
conc_flow  = (s1_tgt ./ total_flow) * 2.0;

valid_c   = conc_flow >= 0 & conc_flow <= 2.5;
flow_t    = flow_t(valid_c);
conc_flow = conc_flow(valid_c);

%% =========================
% 3. EN iYi OFFSET BUL
%    pp.mat 0-88 dk, akis verisi 0-421 dk
%    Sweep zamanlarini farkli offset'lerle kaydir,
%    en iyi R2'yi veren offset'i sec
% =========================
fprintf('Offset optimizasyonu basliyor...\n');

% Test edilecek offset araligI: 0 ile (flow_max - sweep_max) saniye arasi
% Her adim: 5 dakika (300 saniye)
offset_max  = max(flow_t) - max(sweep_time_s);
offset_step = 300;   % 5 dakika adimlarla
offsets     = 0 : offset_step : floor(offset_max/offset_step)*offset_step;

best_R2     = -inf;
best_offset = 0;
R2_per_offset = zeros(size(offsets));

% Ozellikleri hazirla
sigma_X = std(X);
X(:, sigma_X == 0) = [];
mu_X    = mean(X);
sigma_X = std(X);
sigma_X(sigma_X == 0) = 1;
Xz = (X - mu_X) ./ sigma_X;

n = size(Xz, 1);

for oi = 1:length(offsets)
    offset = offsets(oi);

    % Sweep zamanlarini offset ile kaydir
    shifted_time = sweep_time_s + offset;

    % Her sweep icin konsantrasyon interpolasyonu
    y_try = interp1(flow_t, conc_flow, shifted_time, 'linear', NaN);

    % NaN olmayan sweepleri al
    valid_sw = ~isnan(y_try);
    if sum(valid_sw) < 20, continue; end

    Xv = Xz(valid_sw, :);
    yv = y_try(valid_sw);

    % Konsantrasyonda yeterli cesitlilik var mi?
    if std(yv) < 0.05, continue; end

    % Basit PLS (3 bilesken) ile hizli R2 hesapla
    mu_yv = mean(yv);
    yv_c  = yv - mu_yv;

    rng(42);
    nv     = length(yv);
    idx_r  = randperm(nv);
    nTrain = round(0.7 * nv);
    tr = idx_r(1:nTrain);
    te = idx_r(nTrain+1:end);

    try
        [B_tmp, ~, ~] = nipals_pls(Xv(tr,:), yv_c(tr), 3);
        y_pred_tmp = Xv(te,:) * B_tmp + mu_yv;
        ss_res = sum((yv(te) - y_pred_tmp).^2);
        ss_tot = sum((yv(te) - mean(yv(te))).^2);
        R2_tmp = 1 - ss_res / ss_tot;
    catch
        R2_tmp = -inf;
    end

    R2_per_offset(oi) = R2_tmp;

    if R2_tmp > best_R2
        best_R2     = R2_tmp;
        best_offset = offset;
    end
end

fprintf('En iyi offset: %.0f saniye (%.1f dakika)\n', best_offset, best_offset/60);
fprintf('En iyi R2    : %.4f\n', best_R2);

%% =========================
% 4. EN iYi OFFSET iLE FINAL y VEKTORU
% =========================
shifted_time_final = sweep_time_s + best_offset;
y = interp1(flow_t, conc_flow, shifted_time_final, 'linear', 'extrap');
y = max(0, min(2.0, y));

fprintf('\nFinal konsantrasyon araligi: %.3f%% - %.3f%%\n', min(y), max(y));
fprintf('Konsantrasyon std: %.4f\n', std(y));

%% =========================
% 5. EGITIM / TEST BOLUMU
% =========================
n   = length(y);
mu_y = mean(y);
y_c  = y - mu_y;

rng(42);
idx_r  = randperm(n);
nTrain = round(0.7 * n);
tr_idx = idx_r(1:nTrain);
te_idx = idx_r(nTrain+1:end);

Xtrain = Xz(tr_idx,:);
ytrain = y_c(tr_idx);
Xtest  = Xz(te_idx,:);
ytest  = y(te_idx);

%% =========================
% 6. PCA (manuel)
% =========================
C = (Xtrain' * Xtrain) / (size(Xtrain,1) - 1);
[coeff, eigVal_mat] = eig(C);
latent = diag(eigVal_mat);
[latent, sidx] = sort(latent, 'descend');
coeff = coeff(:, sidx);

score_tr      = Xtrain * coeff;
explained_cum = cumsum(latent / sum(latent));
k = find(explained_cum >= 0.95, 1);
if isempty(k), k = size(coeff,2); end
fprintf('  PCA bilesken sayisi: %d\n', k);

Xtrain_pca = score_tr(:, 1:k);
Xtest_pca  = (Xtest - mean(Xtrain)) * coeff(:, 1:k);

%% =========================
% 7. FINAL PLS
% =========================
nComp = min(k, 10);
[B_pls, ~, pctvar] = nipals_pls(Xtrain_pca, ytrain, nComp);
y_pred_pls = Xtest_pca * B_pls + mu_y;

R2_pls   = 1 - sum((ytest - y_pred_pls).^2) / sum((ytest - mean(ytest)).^2);
RMSE_pls = sqrt(mean((ytest - y_pred_pls).^2));
LOD_pls  = 3 * RMSE_pls;

fprintf('\n=========================\n');
fprintf('FINAL PLS (offset=%.0fs)\n', best_offset);
fprintf('R2   = %.4f\n',   R2_pls);
fprintf('RMSE = %.4f %%\n', RMSE_pls);
fprintf('LOD  = %.4f %%\n', LOD_pls);
fprintf('=========================\n');

%% =========================
% 8. GRAFIKLER
% =========================

% --- R2 vs Offset ---
figure('Name','R2 vs Offset','Position',[100 100 600 400]);
plot(offsets/60, R2_per_offset, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor','b');
hold on;
plot(best_offset/60, best_R2, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Offset (dakika)');
ylabel('R^2');
title('PLS R^2 — Offset Optimizasyonu (data4)');
legend('R^2', sprintf('En iyi: %.0f dk', best_offset/60));
grid on;

% --- PLS Parity Plot ---
figure('Name','PLS Sonuclari','Position',[100 550 1200 400]);

subplot(1,3,1);
plot(ytest, 'b', 'LineWidth', 1.5); hold on;
plot(y_pred_pls, 'r--', 'LineWidth', 1.5);
legend('Gercek','Tahmin');
title(sprintf('PLS Zaman Serisi (R^2=%.4f)', R2_pls));
xlabel('Ornek'); ylabel('Konsantrasyon (%)');
grid on;

subplot(1,3,2);
scatter(ytest, y_pred_pls, 40, 'filled', 'MarkerFaceColor',[0.2 0.4 0.8]);
hold on;
ref = linspace(min(ytest), max(ytest), 100);
plot(ref, ref, 'r--', 'LineWidth', 1.5);
xlabel('Gercek (%)'); ylabel('Tahmin (%)');
title(sprintf('PLS Parity Plot\nR^2=%.4f, LOD=%.4f%%', R2_pls, LOD_pls));
grid on;

subplot(1,3,3);
plot(ytest - y_pred_pls, 'k');
yline(0, 'r--');
title('PLS Kalintilari');
xlabel('Ornek'); ylabel('Hata (%)');
grid on;

sgtitle(sprintf('LAG-Optimize PLS — data4 (Offset=%.0f dk)', best_offset/60));

% --- Konsantrasyon Zaman Serisi ---
figure('Name','Konsantrasyon Hizalama','Position',[100 100 900 380]);
plot(sweep_time_s/60, y, 'r-', 'LineWidth', 1.5);
xlabel('Sweep Zamani (dakika)');
ylabel('Konsantrasyon (%)');
title(sprintf('Hizalanmis Konsantrasyon (Offset=%.0f dk)', best_offset/60));
grid on;


%% =========================================================
%  YARDIMCI FONKSiYON: NIPALS-PLS
% =========================================================
function [B, T, pctvar] = nipals_pls(X, y, nComp)
    [n, p]  = size(X);
    X0 = X;  y0 = y;
    y_ss = sum(y.^2);

    W = zeros(p, nComp);
    P = zeros(p, nComp);
    Q = zeros(nComp, 1);
    T = zeros(n, nComp);
    pctvar = zeros(1, nComp);

    for a = 1:nComp
        w = X0' * y0;
        w = w / norm(w);
        t = X0 * w;
        t_norm = t' * t;
        p_vec  = X0' * t / t_norm;
        q      = (y0' * t)  / t_norm;

        X0 = X0 - t * p_vec';
        y0 = y0 - t * q;

        W(:,a) = w;  P(:,a) = p_vec;
        Q(a)   = q;  T(:,a) = t;
        pctvar(a) = 100 * sum((t*q).^2) / y_ss;
    end

    B = W * ((P'*W) \ Q);
end