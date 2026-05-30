clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

load(fullfile(scriptDir, 'pp.mat'));

figDir = fullfile(scriptDir, 'figures');
resultDir = fullfile(scriptDir, 'results');

if ~exist(figDir, 'dir')
    mkdir(figDir);
end
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

sweep = (1:numel(R_peak))';

features = [R_peak(:), X1_peak(:), X2_peak(:), Z_abs_peak(:), ...
            Phase_peak(:), Y_abs_peak(:), G_peak(:), B_peak(:)];

featureNames = {'R_peak','X1_peak','X2_peak','Z_abs_peak', ...
                'Phase_peak','Y_abs_peak','G_peak','B_peak'};

% Exclude pump degradation after sweep 1785
valid = sweep < 1785;

featuresValid = features(valid, :);
sweepValid = sweep(valid);

% Segment valid sweeps into 21 concentration plateaus: 0.0 to 2.0%
nLevels = 21;
concLevels = (0:0.1:2.0)';
edges = round(linspace(min(sweepValid), max(sweepValid) + 1, nLevels + 1));

plateauID = nan(size(sweep));
concSweep = nan(size(sweep));

for k = 1:nLevels
    idx = sweep >= edges(k) & sweep < edges(k+1) & valid;
    plateauID(idx) = k;
    concSweep(idx) = concLevels(k);
end

% Figure 1: time/sweep series
figure('Color','w');
plot(sweep, R_peak, 'LineWidth', 1.2);
hold on;
xline(1785, '--r', 'Pump failure onset');
xlabel('Sweep number');
ylabel('R peak frequency (Hz)');
title('R peak resonance frequency over experiment');
grid on;
saveas(gcf, fullfile(figDir, 'r_peak_timeseries.png'));

% Plateau statistics
meanVals = nan(nLevels, numel(featureNames));
stdVals = nan(nLevels, numel(featureNames));

for k = 1:nLevels
    idx = plateauID == k;
    meanVals(k, :) = mean(features(idx, :), 1, 'omitnan');
    stdVals(k, :) = std(features(idx, :), 0, 1, 'omitnan');
end

meanTable = array2table(meanVals, 'VariableNames', featureNames);
stdTable = array2table(stdVals, 'VariableNames', strcat(featureNames, '_std'));

plateauStats = [table(concLevels, 'VariableNames', {'Concentration_percent'}), ...
                meanTable, stdTable];

writetable(plateauStats, fullfile(resultDir, 'plateau_statistics.csv'));

% Calibration analysis for three representative features
selectedIdx = [1 2 7]; % R_peak, X1_peak, G_peak
calibResults = table();

figure('Color','w');
tiledlayout(1, 3, 'TileSpacing', 'compact');

for j = 1:numel(selectedIdx)
    col = selectedIdx(j);

    x = concLevels;
    y = meanVals(:, col);
    yerr = stdVals(:, col);

    [intercept, sensitivity, R2, xfit, yfit, yci] = simpleLinearFit(x, y);

    noiseFloor = yerr(1);
    LOD = 3 * noiseFloor / abs(sensitivity);

    calibResults.Feature(j, 1) = string(featureNames{col});
    calibResults.Sensitivity_Hz_per_percent(j, 1) = sensitivity;
    calibResults.Intercept_Hz(j, 1) = intercept;
    calibResults.R2(j, 1) = R2;
    calibResults.NoiseFloor_Hz(j, 1) = noiseFloor;
    calibResults.LOD_percent(j, 1) = LOD;

    nexttile;
    errorbar(x, y, yerr, 'o', 'LineWidth', 1.1);
    hold on;

    plot(xfit, yfit, 'r-', 'LineWidth', 1.5);
    plot(xfit, yci, 'r--', 'LineWidth', 0.8);

    xlabel('Glycerol concentration (% w/v)');
    ylabel('Frequency (Hz)');
    title(strrep(featureNames{col}, '_', '\_'));
    grid on;
end

saveas(gcf, fullfile(figDir, 'calibration_curves.png'));
writetable(calibResults, fullfile(resultDir, 'calibration_results.csv'));

% PCA on plateau mean feature matrix using SVD, avoiding toolbox dependency
X = standardizeColumns(meanVals);
[U, S, coeff] = svd(X, 'econ');
score = U * S;
latent = diag(S).^2 ./ (size(X, 1) - 1);
explained = 100 * latent ./ sum(latent);

pcaTable = table((1:numel(explained))', explained, cumsum(explained), ...
    'VariableNames', {'PC', 'Explained_percent', 'Cumulative_percent'});
writetable(pcaTable, fullfile(resultDir, 'pca_explained_variance.csv'));

figure('Color','w');
scatter(score(:,1), score(:,2), 50, concLevels, 'filled');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA of standardized impedance features');
cb = colorbar;
cb.Label.String = 'Glycerol concentration (% w/v)';
grid on;
saveas(gcf, fullfile(figDir, 'pca_result.png'));

% Simple multivariate regression using first two PCs
Xreg = score(:, 1:2);
yreg = concLevels;

[pcrIntercept, pcrBeta, pcrR2, yhat] = simpleMultipleFit(Xreg, yreg);
rmsePCR = sqrt(mean((yreg - yhat).^2));

figure('Color','w');
plot(yreg, yhat, 'o', 'LineWidth', 1.2);
hold on;
plot([0 2], [0 2], 'k--');
xlabel('Reference concentration (% w/v)');
ylabel('Predicted concentration (% w/v)');
title(sprintf('PCR prediction, RMSE = %.4f %%', rmsePCR));
grid on;
saveas(gcf, fullfile(figDir, 'pcr_prediction.png'));

multivarResults = table(rmsePCR, pcrR2, ...
    'VariableNames', {'PCR_RMSE_percent', 'PCR_R2'});
writetable(multivarResults, fullfile(resultDir, 'multivariate_results.csv'));
disp('Calibration results:');
disp(calibResults);

disp('PCA explained variance:');
disp(pcaTable(1:3, :));

disp('Done. Figures saved in /figures and CSV files saved in /results.');
function [intercept, slope, R2, xfit, yfit, yci] = simpleLinearFit(x, y)
    x = x(:);
    y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    A = [ones(size(x)), x];
    beta = A \ y;
    intercept = beta(1);
    slope = beta(2);

    yhat = A * beta;
    residual = y - yhat;
    sse = sum(residual.^2);
    sst = sum((y - mean(y)).^2);
    R2 = 1 - sse / sst;

    xfit = linspace(min(x), max(x), 100)';
    yfit = intercept + slope .* xfit;

    n = numel(x);
    s2 = sse / max(n - 2, 1);
    xbar = mean(x);
    sxx = sum((x - xbar).^2);
    sePred = sqrt(s2 .* (1 + 1/n + ((xfit - xbar).^2 ./ sxx)));
    yci = [yfit - 2 * sePred, yfit + 2 * sePred];
end

function Xs = standardizeColumns(X)
    mu = mean(X, 1, 'omitnan');
    sigma = std(X, 0, 1, 'omitnan');
    sigma(sigma == 0 | ~isfinite(sigma)) = 1;
    Xs = (X - mu) ./ sigma;
    Xs(~isfinite(Xs)) = 0;
end

function [intercept, slopes, R2, yhat] = simpleMultipleFit(X, y)
    y = y(:);
    keep = all(isfinite(X), 2) & isfinite(y);
    X = X(keep, :);
    y = y(keep);

    A = [ones(size(X, 1), 1), X];
    beta = A \ y;
    intercept = beta(1);
    slopes = beta(2:end);
    yhat = A * beta;

    sse = sum((y - yhat).^2);
    sst = sum((y - mean(y)).^2);
    R2 = 1 - sse / sst;
end
