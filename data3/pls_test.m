clearvars
clc
dataDir = fileparts(mfilename('fullpath'));
dataRoot = fileparts(dataDir);
%% =========================
% 1. LOAD FEATURE DATA
% =========================
if isfile(fullfile(dataRoot,'3rd_fit_v3.mat'))
    load('3rd_fit_v3.mat','T_all');
    X = table2array(T_all);
end

%% =========================
% 2. LOAD TARGET DATA
% =========================
baseDir = fileparts(mfilename('fullpath'));
dataRoot = fileparts(baseDir);
dataFile = fullfile(dataRoot, 'IWIS2025flow', 'data.txt');
opts = detectImportOptions(dataFile);
opts.VariableNamingRule = 'modify';
Traw = readtable(dataFile, opts);
s1 = Traw{:,4};
s2 = Traw{:,5};
conc = (s1 ./ (s1+s2)) * 2.0;

%% =========================
% 3. ALIGN DATA
% =========================
pp = load("pp.mat");
t_sweep_s = pp.R_time(:) * 60;
t_flow = Traw{:,2};
t_flow = t_flow - t_flow(1);
y = interp1(t_flow, conc, t_sweep_s, "linear", "extrap");
y = y(:);

% --- Görsel kontrol ---
figure('Name','Hizalama Kontrolü');
subplot(2,1,1); plot(t_flow, conc); title('Ham Konsantrasyon'); xlabel('Zaman (s)'); ylabel('Conc');
subplot(2,1,2); plot(t_sweep_s, y);  title('İnterpolasyon Sonrası'); xlabel('Zaman (s)'); ylabel('Conc');

n = size(X,1);

%% =========================
% 4. REMOVE CONSTANT / LOW-VARIANCE FEATURES
% =========================
sigma_raw = std(X);
X(:, sigma_raw == 0) = [];

%% =========================
% 5. STANDARDIZATION
% =========================
mu_x    = mean(X);
sigma_x = std(X);
sigma_x(sigma_x == 0) = 1;
Xz = (X - mu_x) ./ sigma_x;

%% =========================
% 6. KORELASYON ANALİZİ (Tanısal)
% =========================
corr_vals = corr(Xz, y);
[sorted_corr, corr_order] = sort(abs(corr_vals), 'descend');
fprintf('=========================\n');
fprintf('En yüksek 10 özellik korelasyonu:\n');
for i = 1:min(10, length(sorted_corr))
    fprintf('  Özellik %d: r = %.4f\n', corr_order(i), sorted_corr(i));
end
fprintf('=========================\n');

%% =========================
% 7. TRAIN / TEST SPLIT
% =========================
rng(1);
idx      = randperm(n);
trainRatio = 0.7;
nTrain   = round(trainRatio * n);
trainIdx = idx(1:nTrain);
testIdx  = idx(nTrain+1:end);

Xtrain = Xz(trainIdx,:);
ytrain = y(trainIdx);
Xtest  = Xz(testIdx,:);
ytest  = y(testIdx);

%% =========================
% 8. PCA — Düzeltilmiş Projeksiyon
% =========================
mu_train    = mean(Xtrain);          % PCA için eğitim ortalaması
Xtrain_c    = Xtrain - mu_train;
Xtest_c     = Xtest  - mu_train;     

[coeff, score, latent] = pca(Xtrain_c);
explained = cumsum(latent / sum(latent));

k95  = find(explained >= 0.95, 1);
k99  = find(explained >= 0.99, 1);
fprintf('PCA: %%95 için %d bileşen, %%99 için %d bileşen\n', k95, k99);

% %99 ile daha fazla bilgi tut
k = k99;
Xtrain_pca = score(:,1:k);
Xtest_pca  = Xtest_c * coeff(:,1:k);   

%% =========================
% 9. MODEL KARŞILAŞTIRMASI
% =========================
results = struct();

% --- A) SVR (Hiperparametre Optimizasyonlu) ---
fprintf('\nSVR eğitimi (optimize ediliyor)...\n');
mdl_svr = fitrsvm(Xtrain_pca, ytrain, ...
    'KernelFunction', 'gaussian', ...
    'OptimizeHyperparameters', {'BoxConstraint','Epsilon','KernelScale'}, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'MaxObjectiveEvaluations', 40, ...
        'Repartition', true, ...
        'ShowPlots', false, ...
        'Verbose', 0), ...
    'Standardize', true);

y_pred_svr = predict(mdl_svr, Xtest_pca);
results.SVR.R2   = 1 - sum((ytest - y_pred_svr).^2) / sum((ytest - mean(ytest)).^2);
results.SVR.RMSE = sqrt(mean((ytest - y_pred_svr).^2));
results.SVR.pred = y_pred_svr;

% --- B) Random Forest (PCA'sız) ---
fprintf('Random Forest eğitimi...\n');
mdl_rf = fitrensemble(Xtrain, ytrain, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 200, ...
    'Learners', templateTree('MinLeafSize', 5));

y_pred_rf = predict(mdl_rf, Xtest);
results.RF.R2   = 1 - sum((ytest - y_pred_rf).^2) / sum((ytest - mean(ytest)).^2);
results.RF.RMSE = sqrt(mean((ytest - y_pred_rf).^2));
results.RF.pred = y_pred_rf;

% --- C) Gradient Boosting ---
fprintf('Gradient Boosting eğitimi...\n');
mdl_gb = fitrensemble(Xtrain, ytrain, ...
    'Method', 'LSBoost', ...
    'NumLearningCycles', 200, ...
    'LearnRate', 0.05, ...
    'Learners', templateTree('MaxNumSplits', 4));

y_pred_gb = predict(mdl_gb, Xtest);
results.GB.R2   = 1 - sum((ytest - y_pred_gb).^2) / sum((ytest - mean(ytest)).^2);
results.GB.RMSE = sqrt(mean((ytest - y_pred_gb).^2));
results.GB.pred = y_pred_gb;

% --- D) GPR (Gaussian Process) —---
fprintf('GPR eğitimi...\n');
mdl_gpr = fitrgp(Xtrain_pca, ytrain, ...
    'KernelFunction', 'squaredexponential', ...
    'Standardize', true, ...
    'OptimizeHyperparameters', 'auto', ...
    'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0));

y_pred_gpr = predict(mdl_gpr, Xtest_pca);
results.GPR.R2   = 1 - sum((ytest - y_pred_gpr).^2) / sum((ytest - mean(ytest)).^2);
results.GPR.RMSE = sqrt(mean((ytest - y_pred_gpr).^2));
results.GPR.pred = y_pred_gpr;

%% =========================
% 10. SONUÇLAR
% =========================
model_names = fieldnames(results);
fprintf('\n=========================================\n');
fprintf('%-10s  %8s  %8s\n', 'MODEL', 'R²', 'RMSE');
fprintf('-----------------------------------------\n');
for i = 1:numel(model_names)
    nm = model_names{i};
    fprintf('%-10s  %8.4f  %8.4f\n', nm, results.(nm).R2, results.(nm).RMSE);
end
fprintf('=========================================\n');

% En iyi modeli bul
r2_vals = cellfun(@(nm) results.(nm).R2, model_names);
[~, best_idx] = max(r2_vals);
best_name = model_names{best_idx};
best_pred = results.(best_name).pred;
fprintf('En iyi model: %s (R² = %.4f)\n\n', best_name, results.(best_name).R2);

%% =========================
% 11. PLOTS — En İyi Model
% =========================
% Zaman serisi karşılaştırması
figure('Name', sprintf('En İyi Model: %s', best_name));
plot(ytest, 'b', 'LineWidth', 1.5); hold on;
plot(best_pred, 'r--', 'LineWidth', 1.5);
legend('Gerçek', 'Tahmin');
title(sprintf('%s Model Performansı (R²=%.4f)', best_name, results.(best_name).R2));
xlabel('Örnek'); ylabel('Akış Hızı'); grid on;

% Parity plot
figure('Name','Parity Plot');
scatter(ytest, best_pred, 40, 'filled', 'MarkerFaceAlpha', 0.7);
hold on; refline(1,0);
xlabel('Gerçek Değerler'); ylabel('Tahmin Değerleri');
title(sprintf('Parity Plot — %s', best_name)); grid on;

% Tüm modeller karşılaştırma
figure('Name','Model Karşılaştırması');
colors = lines(numel(model_names));
plot(ytest, 'k-', 'LineWidth', 2, 'DisplayName', 'Gerçek'); hold on;
for i = 1:numel(model_names)
    nm = model_names{i};
    plot(results.(nm).pred, '--', 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s (R²=%.3f)', nm, results.(nm).R2));
end
legend('Location','best'); title('Model Karşılaştırması');
xlabel('Örnek'); ylabel('Akış Hızı'); grid on;

% Residual analizi — en iyi model
figure('Name','Residual Analizi');
residuals = ytest - best_pred;
subplot(1,2,1);
plot(residuals, 'Color', [0.2 0.5 0.8]);
title('Residuals'); xlabel('Örnek'); ylabel('Hata'); grid on;
yline(0, 'r--');
subplot(1,2,2);
histogram(residuals, 20, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'white');
title('Residual Dağılımı'); xlabel('Hata'); ylabel('Frekans'); grid on;

% PCA Açıklanan Varyans
figure('Name','PCA Açıklanan Varyans');
plot(explained * 100, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xline(k95, 'r--', '95%'); xline(k99, 'g--', '99%');
xlabel('Bileşen Sayısı'); ylabel('Kümülatif Açıklanan Varyans (%)');
title('PCA Scree Plot'); grid on; ylim([0 101]);