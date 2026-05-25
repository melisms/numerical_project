clearvars
clc

%% =========================
% 1. LOAD FEATURE DATA
% =========================
load('3rd_fit.mat','T');
X = table2array(T);

%% =========================
% 2. LOAD TARGET DATA
% =========================
opts = detectImportOptions('/Users/melisamuslu/Documents/MATLAB/data/IWIS2025flow/data.txt');
opts.VariableNamingRule = 'modify';
Traw = readtable('/Users/melisamuslu/Documents/MATLAB/data/IWIS2025flow/data.txt', opts);

y = Traw{:,4};   % s1(Read)

%% =========================
% 3. ALIGN DATA (CRITICAL)
% =========================
n = min(size(X,1), length(y));

X = X(1:n,:);
y = y(1:n,:);

%% =========================
% 4. REMOVE CONSTANT FEATURES
% =========================
sigma = std(X);
X(:, sigma == 0) = [];

%% =========================
% 5. STANDARDIZATION
% =========================
mu = mean(X);
sigma = std(X);
sigma(sigma==0) = 1;

Xz = (X - mu) ./ sigma;

%% =========================
% 6. TRAIN / TEST SPLIT
% =========================
rng(1); % reproducibility
idx = randperm(n);

trainRatio = 0.7;
nTrain = round(trainRatio * n);

trainIdx = idx(1:nTrain);
testIdx  = idx(nTrain+1:end);

Xtrain = Xz(trainIdx,:);
ytrain = y(trainIdx);

Xtest  = Xz(testIdx,:);
ytest  = y(testIdx);

%% =========================
% 7. OPTIONAL: PCA (feature compression)
% =========================
[coeff, score, latent] = pca(Xtrain);

explained = cumsum(latent / sum(latent));
k = find(explained >= 0.95, 1);

Xtrain_pca = score(:,1:k);
Xtest_pca = (Xtest - mean(Xtrain)) * coeff(:,1:k);

%% =========================
% 8. NONLINEAR MODEL (SVR)
% =========================
mdl = fitrsvm(Xtrain_pca, ytrain, ...
    'KernelFunction','gaussian', ...
    'Standardize',true);

%% =========================
% 9. PREDICTION
% =========================
y_pred = predict(mdl, Xtest_pca);

%% =========================
% 10. PERFORMANCE METRICS
% =========================
R2 = 1 - sum((ytest - y_pred).^2) / sum((ytest - mean(ytest)).^2);
RMSE = sqrt(mean((ytest - y_pred).^2));

fprintf('=========================\n');
fprintf('MODEL PERFORMANCE\n');
fprintf('R2   = %.4f\n', R2);
fprintf('RMSE = %.4f\n', RMSE);
fprintf('=========================\n');

%% =========================
% 11. PLOTS
% =========================

% Time series comparison
figure;
plot(ytest,'b','LineWidth',1.5); hold on;
plot(y_pred,'r--','LineWidth',1.5);
legend('True','Predicted');
title('SVR Model Performance');
xlabel('Samples');
ylabel('Flow Rate');
grid on;

% Parity plot
figure;
scatter(ytest, y_pred, 25, 'filled');
xlabel('True Values');
ylabel('Predicted Values');
title('Parity Plot');
grid on;
refline(1,0);

%% =========================
% 12. RESIDUAL ANALYSIS
% =========================
figure;
residuals = ytest - y_pred;
plot(residuals);
title('Residuals');
xlabel('Samples');
ylabel('Error');
grid on;