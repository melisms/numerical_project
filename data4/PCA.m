%% =========================================================
%  PCA_data4.m
%  3rd_fit_data4.mat uzerinde PCA analizi
%
%  On kosul: mdpi_poly_data4.m calistirilmis olmali
%            (3rd_fit_data4.mat olusturulmus olmali)
% =========================================================
clearvars
clc

%% =========================
% 1. VERI YUKLE
% =========================
load('3rd_fit_data4.mat', 'T');
X = table2array(T);
featureNames = T.Properties.VariableNames;

fprintf('Yuklendi: %d sweep x %d ozellik\n', size(X,1), size(X,2));

%% =========================
% 2. STANDARDIZE ET
% =========================
mu    = mean(X, 1);
sigma = std(X, 0, 1);
sigma(sigma == 0) = 1;
Xz = (X - mu) ./ sigma;

%% =========================
% 3. KOVARYANS MATRiSi
% =========================
C = (Xz' * Xz) / (size(Xz,1) - 1);

%% =========================
% 4. OZVEKTOR / OZDEGER
% =========================
[eigVec, eigVal_mat] = eig(C);
eigVal = diag(eigVal_mat);

[sortedEigVal, idx] = sort(eigVal, 'descend');
eigVec = eigVec(:, idx);

%% =========================
% 5. PCA SKORLARI
% =========================
score = Xz * eigVec;

%% =========================
% 6. ACIKLANAN VARYANS
% =========================
explained  = 100 * sortedEigVal / sum(sortedEigVal);
cumulative = cumsum(explained);

fprintf('\n  PC1 aciklanan varyans : %.2f%%\n', explained(1));
fprintf('  PC2 aciklanan varyans : %.2f%%\n', explained(2));
fprintf('  %%95 icin gereken PC  : %d\n', find(cumulative >= 95, 1));

%% =========================
% 7. GRAFIKLER
% =========================

% --- PCA Score Plot ---
figure('Name','PCA Score Plot','Position',[100 100 550 450]);
scatter(score(:,1), score(:,2), 40, 'filled');
xlabel('PC1'); ylabel('PC2');
title('PCA Score Plot — data4');
grid on;

% --- Scree Plot ---
figure('Name','Scree Plot','Position',[670 100 550 450]);
bar(explained(1:min(10,end)), 'FaceColor', [0.2 0.5 0.8]);
xlabel('Temel Bilesken');
ylabel('Aciklanan Varyans (%)');
title('Scree Plot — data4');
grid on;

% --- Kumulatif Varyans ---
figure('Name','Kumulatif Varyans','Position',[100 580 550 400]);
plot(cumulative, '-o', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
xlabel('Bilesken Sayisi');
ylabel('Kumulatif Varyans (%)');
title('Kumulatif Aciklanan Varyans — data4');
yline(95, '--k', '%95', 'LabelHorizontalAlignment','left');
grid on;

% --- PC1 Loadings ---
figure('Name','PC1 Loadings','Position',[670 580 700 400]);
bar(eigVec(:,1), 'FaceColor', [0.4 0.7 0.3]);
title('PC1 Yukleme Grafigi — data4');
xticks(1:length(featureNames));
xticklabels(featureNames);
xtickangle(45);
ylabel('Yukleme Degeri');
grid on;

% --- PC2 Loadings ---
figure('Name','PC2 Loadings','Position',[100 580 700 400]);
bar(eigVec(:,2), 'FaceColor', [0.9 0.5 0.2]);
title('PC2 Yukleme Grafigi — data4');
xticks(1:length(featureNames));
xticklabels(featureNames);
xtickangle(45);
ylabel('Yukleme Degeri');
grid on;

%% =========================
% 8. EN ONEMLI OZELLIKLER
% =========================
[~, idx1] = sort(abs(eigVec(:,1)), 'descend');
[~, idx2] = sort(abs(eigVec(:,2)), 'descend');

disp('===== EN ONEMLI OZELLIKLER (PC1) =====');
for i = 1:min(10, length(idx1))
    fprintf('%s : %.4f\n', featureNames{idx1(i)}, eigVec(idx1(i),1));
end

disp('===== EN ONEMLI OZELLIKLER (PC2) =====');
for i = 1:min(10, length(idx2))
    fprintf('%s : %.4f\n', featureNames{idx2(i)}, eigVec(idx2(i),2));
end

% PCA sonuclarini kaydet (pls_data4.m icin)
save('pca_result_data4.mat', 'eigVec', 'score', 'explained', 'cumulative', 'mu', 'sigma');
fprintf('\npca_result_data4.mat kaydedildi.\n');