clearvars
clc

%% =========================
% 1. LOAD DATA
% =========================
load('3rd_fit.mat','T');

X = table2array(T);

featureNames = T.Properties.VariableNames;

%% =========================
% 2. STANDARDIZATION
% =========================
mu = mean(X,1);
sigma = std(X,0,1);
sigma(sigma == 0) = 1;

Xz = (X - mu) ./ sigma;

%% =========================
% 3. COVARIANCE MATRIX
% =========================
C = (Xz' * Xz) / (size(Xz,1)-1);

%% =========================
% 4. EIGEN DECOMPOSITION (FIXED)
% =========================
[eigVec, eigVal] = eig(C);

eigVal = diag(eigVal);

% sort descending
[sortedEigVal, idx] = sort(eigVal, 'descend');
eigVec = eigVec(:, idx);

%% =========================
% 5. PCA SCORES
% =========================
score = Xz * eigVec;

%% =========================
% 6. EXPLAINED VARIANCE
% =========================
explained = 100 * sortedEigVal / sum(sortedEigVal);
cumulative = cumsum(explained);

%% =========================
% 7. PLOTS
% =========================

% --- PCA scatter
figure;
scatter(score(:,1), score(:,2), 40, 'filled')
xlabel('PC1')
ylabel('PC2')
title('PCA Score Plot')
grid on

% --- Scree plot
figure;
bar(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('Scree Plot')
grid on

% --- Cumulative variance
figure;
plot(cumulative, '-o','LineWidth',2)
xlabel('Number of Components')
ylabel('Cumulative Variance (%)')
title('Cumulative Explained Variance')
grid on

%% =========================
% 8. FEATURE IMPORTANCE (LOADINGS)
% =========================

PC1_loading = abs(eigVec(:,1));
PC2_loading = abs(eigVec(:,2));

[sorted_PC1, idx1] = sort(PC1_loading, 'descend');
[sorted_PC2, idx2] = sort(PC2_loading, 'descend');

%% =========================
% 9. DISPLAY TOP FEATURES
% =========================

disp('===== TOP FEATURES (PC1) =====')
for i = 1:min(10,length(idx1))
    fprintf('%s : %.4f\n', featureNames{idx1(i)}, sorted_PC1(i));
end

disp('===== TOP FEATURES (PC2) =====')
for i = 1:min(10,length(idx2))
    fprintf('%s : %.4f\n', featureNames{idx2(i)}, sorted_PC2(i));
end

%% =========================
% 10. LOADINGS VISUALIZATION
% =========================

figure;
bar(eigVec(:,1))
title('PC1 Loadings')
xticks(1:length(featureNames))
xticklabels(featureNames)
xtickangle(45)

figure;
bar(eigVec(:,2))
title('PC2 Loadings')
xticks(1:length(featureNames))
xticklabels(featureNames)
xtickangle(45)