%% =========================================================
%  data3_pcr_clean.m
%  data3 Multivariate Analizi: PCA + PCR + VIF
%  (Kalibrasyon CIKARILDI - sweep/akis zaman uyusmazligi nedeniyle)
%
%  KULLANIM:
%  - Bu dosyayi data3 klasorune koy (pp.mat ile ayni yere)
%  - data.txt, ust klasordeki IWIS2025flow icinde olmali
%  - MATLAB'da Current Folder'i data3 yap, calistir
% =========================================================
clear; clc; close all;
warning('off', 'MATLAB:table:RowsAddedExistingVars');

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

load(fullfile(scriptDir, 'pp.mat'));

resultDir = fullfile(scriptDir, 'results');
figDir    = fullfile(scriptDir, 'figures');
if ~exist(resultDir, 'dir'), mkdir(resultDir); end
if ~exist(figDir,    'dir'), mkdir(figDir);    end

% =====================
% 1. EMPEDANS OZELLIKLERI - pp.mat
% =====================
sweep = (1:numel(R_peak))';

features = [R_peak(:), X1_peak(:), X2_peak(:), Z_abs_peak(:), ...
            Phase_peak(:), Y_abs_peak(:), G_peak(:), B_peak(:)];

featureNames = {'R_peak','X1_peak','X2_peak','Z_abs_peak', ...
                'Phase_peak','Y_abs_peak','G_peak','B_peak'};

% Pompa arizasi - sweep 1785 sonrasi haric
valid      = sweep < 1785;
sweepValid = sweep(valid);

% =====================
% 2. PLATEAU SEGMENTASYONU (21 seviye)
% =====================
nLevels = 21;
edges   = round(linspace(min(sweepValid), max(sweepValid)+1, nLevels+1));

plateauID = nan(size(sweep));
for k = 1:nLevels
    idx = sweep >= edges(k) & sweep < edges(k+1) & valid;
    plateauID(idx) = k;
end

plateauConc = (0:0.1:2.0)';

% Her plato icin ortalama (gurultu burada yok olur)
meanVals = nan(nLevels, numel(featureNames));
for k = 1:nLevels
    idx = plateauID == k;
    meanVals(k,:) = mean(features(idx,:), 1, 'omitnan');
end

% =====================
% 3. PCA (SVD ile, toolbox yok)
% =====================
X_std = standardizeColumns(meanVals);
[U, S, coeff] = svd(X_std, 'econ');
score     = U * S;
latent    = diag(S).^2 ./ (size(X_std,1)-1);
explained = 100 * latent ./ sum(latent);

pcaTable = table((1:numel(explained))', explained, cumsum(explained), ...
    'VariableNames',{'PC','Explained_percent','Cumulative_percent'});
writetable(pcaTable, fullfile(resultDir,'pca_explained_variance.csv'));

figure('Color','w');
scatter(score(:,1), score(:,2), 50, plateauConc, 'filled');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA of standardized impedance features');
cb = colorbar; cb.Label.String = 'Glycerol concentration (% w/v)';
grid on;
saveas(gcf, fullfile(figDir,'pca_result.png'));

fprintf('PCA explained variance:\n');
disp(pcaTable(1:4,:));

% =====================
% 4. PCR (konsantrasyon tahmini)
% =====================
Xreg = score(:,1:2);
yreg = plateauConc;

[~, ~, pcrR2, yhat] = simpleMultipleFit(Xreg, yreg);
rmsePCR = sqrt(mean((yreg - yhat).^2));

figure('Color','w');
plot(yreg, yhat, 'o', 'LineWidth', 1.2); hold on;
plot([min(plateauConc) max(plateauConc)], ...
     [min(plateauConc) max(plateauConc)], 'k--');
xlabel('Reference concentration (% w/v)');
ylabel('Predicted concentration (% w/v)');
title(sprintf('PCR prediction  |  R^2 = %.4f, RMSE = %.4f %%', pcrR2, rmsePCR));
grid on;
saveas(gcf, fullfile(figDir,'pcr_prediction.png'));

multivarResults = table(rmsePCR, pcrR2, ...
    'VariableNames',{'PCR_RMSE_percent','PCR_R2'});
writetable(multivarResults, fullfile(resultDir,'multivariate_results.csv'));

fprintf('\nPCR R2 = %.4f, RMSE = %.4f %%\n', pcrR2, rmsePCR);

% =====================
% 5. VIF (collinearity, toolbox yok)
% =====================
fprintf('\nVIF Analizi:\n');
fprintf('%-15s %12s\n','Feature','VIF');
fprintf('%s\n', repmat('-',1,28));

vifVals = zeros(1, numel(featureNames));
for f = 1:numel(featureNames)
    y_vif  = X_std(:,f);
    X_vif  = X_std(:, setdiff(1:numel(featureNames), f));
    b_vif  = X_vif \ y_vif;
    yh_vif = X_vif * b_vif;
    ss_res = sum((y_vif - yh_vif).^2);
    ss_tot = sum((y_vif - mean(y_vif)).^2);
    R2_vif = 1 - ss_res / ss_tot;
    vifVals(f) = 1 / max(1 - R2_vif, 1e-10);
    fprintf('%-15s %12.2f\n', featureNames{f}, vifVals(f));
end

vifTable = table(featureNames', vifVals', 'VariableNames',{'Feature','VIF'});
writetable(vifTable, fullfile(resultDir,'vif_results.csv'));

% =====================
% 6. LOD (PCR bazli)
% =====================
LOD_PCR = 3 * rmsePCR;
fprintf('\nPCR LOD = %.4f %%\n', LOD_PCR);

lodTable = table({'PCR (multivariate)'}, LOD_PCR, ...
    'VariableNames',{'Method','LOD_percent'});
writetable(lodTable, fullfile(resultDir,'lod_comparison.csv'));

disp('Done. Figures -> /figures   |   CSV -> /results');

% ============================================================
% LOCAL FUNCTIONS
% ============================================================
function Xs = standardizeColumns(X)
    mu    = mean(X, 1, 'omitnan');
    sigma = std(X,  0, 1, 'omitnan');
    sigma(sigma == 0 | ~isfinite(sigma)) = 1;
    Xs = (X - mu) ./ sigma;
    Xs(~isfinite(Xs)) = 0;
end

function [intercept, slopes, R2, yhat] = simpleMultipleFit(X, y)
    y    = y(:);
    keep = all(isfinite(X),2) & isfinite(y);
    X    = X(keep,:); y = y(keep);
    A    = [ones(size(X,1),1), X];
    beta = A \ y;
    intercept = beta(1);
    slopes    = beta(2:end);
    yhat      = A * beta;
    sse = sum((y - yhat).^2);
    sst = sum((y - mean(y)).^2);
    R2  = 1 - sse/sst;
end