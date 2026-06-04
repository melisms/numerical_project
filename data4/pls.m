%% =========================================================
%  pls_data4.m
%  PLS Regresyonu — data4
%  Hedef degisken: gliserol konsantrasyonu (%)
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
% 2. KONSANTRASYON VEKTORU OLUSTUR
%    pp.mat zaman ekseni (dakika) -> Auto_saved_file4'ten
%    her sweep icin konsantrasyon interpolasyonu
% =========================
pp = load('pp.mat');
sweep_time_min = pp.R_time(:);        % dakika, 698 eleman
sweep_time_s   = sweep_time_min * 60; % saniyeye cevir

flow_path = 'C:\Users\asli\Desktop\numerical_project-main\IWIS2025flow\Auto_saved_file4.txt';
opts = detectImportOptions(flow_path, 'Delimiter', '\t');
opts.VariableNamingRule = 'modify';
opts.Encoding = 'ISO-8859-1';
Traw  = readtable(flow_path, opts);

flow_t  = Traw{:,2};   % zaman (saniye)
s1_tgt  = Traw{:,3};   % s1 hedef akis
s2_tgt  = Traw{:,5};   % s2 hedef akis

% Gecersiz satirlari temizle
ok = ~isnan(flow_t) & ~isnan(s1_tgt) & ~isnan(s2_tgt);
flow_t = flow_t(ok);
s1_tgt = s1_tgt(ok);
s2_tgt = s2_tgt(ok);

% Konsantrasyon hesapla: s1/(s1+s2) * 2.0
total_flow = s1_tgt + s2_tgt;
conc_flow  = (s1_tgt ./ total_flow) * 2.0;

% Saçma degerleri temizle (0-2.5 arasi)
valid_c = conc_flow >= 0 & conc_flow <= 2.5;
flow_t    = flow_t(valid_c);
conc_flow = conc_flow(valid_c);

% Her sweep icin konsantrasyon interpolasyonu
y = interp1(flow_t, conc_flow, sweep_time_s, 'linear', 'extrap');
y = max(0, min(2.0, y));  % 0-2 araligina kist

fprintf('Konsantrasyon araligi: %.3f%% - %.3f%%\n', min(y), max(y));

%% =========================
% 3. UZUNLUK HiZALA
% =========================
n = min(size(X,1), length(y));
X = X(1:n,:);
y = y(1:n);
fprintf('Ornek sayisi: %d\n', n);

%% =========================
% 4. SABiT OZELLIK SIFTIRLA
% =========================
sigma_X = std(X);
X(:, sigma_X == 0) = [];

%% =========================
% 5. STANDARDiZASYON
% =========================
mu_X    = mean(X);
sigma_X = std(X);
sigma_X(sigma_X == 0) = 1;
Xz = (X - mu_X) ./ sigma_X;

mu_y = mean(y);
y_c  = y - mu_y;   % merkeze alinmis y

%% =========================
% 6. EGITIM / TEST BOLUMU
% =========================
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
% 7. PCA ile BOYUT AZALTMA (manuel)
% =========================
C = (Xtrain' * Xtrain) / (size(Xtrain,1) - 1);
[coeff, eigVal_mat] = eig(C);
latent = diag(eigVal_mat);
[latent, sort_idx] = sort(latent, 'descend');
coeff = coeff(:, sort_idx);

score_tr = Xtrain * coeff;
explained_cum = cumsum(latent / sum(latent));
k = find(explained_cum >= 0.95, 1);
if isempty(k), k = size(coeff,2); end
fprintf('  PCA bilesken sayisi: %d (%%95 varyans)\n', k);

Xtrain_pca = score_tr(:, 1:k);
Xtest_pca  = (Xtest - mean(Xtrain)) * coeff(:, 1:k);

%% =========================
% 8. PLS — NIPALS (manuel)
% =========================
nComp = min(k, 10);
[B_pls, ~, pctvar] = nipals_pls(Xtrain_pca, ytrain, nComp);

y_pred_pls = Xtest_pca * B_pls + mu_y;

R2_pls   = 1 - sum((ytest - y_pred_pls).^2) / sum((ytest - mean(ytest)).^2);
RMSE_pls = sqrt(mean((ytest - y_pred_pls).^2));
LOD_pls  = 3 * RMSE_pls;

fprintf('\n=========================\n');
fprintf('PLS MODEL PERFORMANSI\n');
fprintf('R2   = %.4f\n', R2_pls);
fprintf('RMSE = %.4f %%\n', RMSE_pls);
fprintf('LOD  = %.4f %%\n', LOD_pls);
fprintf('=========================\n');

%% =========================
% 9. RIDGE REGRESYON (manuel)
% =========================
lambda = 1e-3;
I_r    = eye(size(Xtrain_pca,2));
B_ridge = (Xtrain_pca' * Xtrain_pca + lambda * I_r) \ (Xtrain_pca' * ytrain);
intercept_r = mu_y - mean(Xtrain_pca) * B_ridge;

y_pred_ridge = Xtest_pca * B_ridge + intercept_r;

R2_ridge   = 1 - sum((ytest - y_pred_ridge).^2) / sum((ytest - mean(ytest)).^2);
RMSE_ridge = sqrt(mean((ytest - y_pred_ridge).^2));

fprintf('\n=========================\n');
fprintf('RIDGE REGRESYON PERFORMANSI\n');
fprintf('R2   = %.4f\n', R2_ridge);
fprintf('RMSE = %.4f %%\n', RMSE_ridge);
fprintf('=========================\n');

%% =========================
% 10. GRAFIKLER — PLS
% =========================
figure('Name','PLS Sonuclari','Position',[100 100 1200 400]);

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
title(sprintf('PLS Parity Plot\nR^2=%.4f, RMSE=%.4f', R2_pls, RMSE_pls));
grid on;

subplot(1,3,3);
plot(ytest - y_pred_pls, 'k');
yline(0, 'r--');
title('PLS Kalintilari');
xlabel('Ornek'); ylabel('Hata (%)');
grid on;

sgtitle('PLS Regresyonu — data4 (Hedef: Konsantrasyon)');

%% =========================
% 11. GRAFIKLER — RIDGE
% =========================
figure('Name','Ridge Sonuclari','Position',[100 550 1200 400]);

subplot(1,3,1);
plot(ytest, 'b', 'LineWidth', 1.5); hold on;
plot(y_pred_ridge, 'r--', 'LineWidth', 1.5);
legend('Gercek','Tahmin');
title(sprintf('Ridge Zaman Serisi (R^2=%.4f)', R2_ridge));
xlabel('Ornek'); ylabel('Konsantrasyon (%)');
grid on;

subplot(1,3,2);
scatter(ytest, y_pred_ridge, 40, 'filled', 'MarkerFaceColor',[0.8 0.3 0.2]);
hold on;
plot(ref, ref, 'r--', 'LineWidth', 1.5);
xlabel('Gercek (%)'); ylabel('Tahmin (%)');
title(sprintf('Ridge Parity Plot\nR^2=%.4f, RMSE=%.4f', R2_ridge, RMSE_ridge));
grid on;

subplot(1,3,3);
plot(ytest - y_pred_ridge, 'k');
yline(0, 'r--');
title('Ridge Kalintilari');
xlabel('Ornek'); ylabel('Hata (%)');
grid on;

sgtitle('Ridge Regresyonu — data4 (Hedef: Konsantrasyon)');

%% =========================
% 12. PLS BILESKEN KATKISI
% =========================
figure('Name','PLS Varyans','Position',[100 100 500 380]);
bar(cumsum(pctvar), 'FaceColor',[0.5 0.7 0.4]);
xlabel('PLS Bilesken');
ylabel('Kumulatif Aciklanan Varyans Y (%)');
title('PLS — Y Varyansi Katkisi');
grid on;


%% =========================================================
%  YARDIMCI FONKSiYON: NIPALS-PLS
% =========================================================
function [B, T, pctvar] = nipals_pls(X, y, nComp)
    [n, p]  = size(X);
    X0 = X;
    y0 = y;
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
        q      = (y0' * t) / t_norm;

        X0 = X0 - t * p_vec';
        y0 = y0 - t * q;

        W(:,a) = w;
        P(:,a) = p_vec;
        Q(a)   = q;
        T(:,a) = t;

        pctvar(a) = 100 * sum((t*q).^2) / y_ss;
    end

    B = W * ((P'*W) \ Q);
end