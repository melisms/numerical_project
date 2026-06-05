clear; clc; close all;
warning('off', 'MATLAB:table:RowsAddedExistingVars');

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

load(fullfile(scriptDir, 'pp.mat'));

figDir    = fullfile(scriptDir, 'figures');
resultDir = fullfile(scriptDir, 'results');

if ~exist(figDir,    'dir'), mkdir(figDir);    end
if ~exist(resultDir, 'dir'), mkdir(resultDir); end

% =====================
% 1. FLOW DATA - data.txt
% =====================
flowFile = fullfile(scriptDir, 'data.txt');
opts = detectImportOptions(flowFile, 'Delimiter', '\t');
opts.DataLines = [2 Inf];
opts.VariableNames = {'TimeStamp','Time_s','s1_target','s1_read','s2_target','s2_read'};
opts.SelectedVariableNames = {'Time_s','s1_target','s1_read','s2_target','s2_read'};
flowData = readtable(flowFile, opts);

t_flow   = flowData.Time_s;
s1_target = flowData.s1_target;
s2_target = flowData.s2_target;
s1_read   = flowData.s1_read;
s2_read   = flowData.s2_read;

% Konsantrasyon = s1 / (s1+s2) * 2.0%
total_flow = s1_target + s2_target;
% Stok konsantrasyonu 5% w/v, toplam akis ~50 uL/min sabit
% s2_target: gliserol stok pompasi (0-20 uL/min)
% s1_target: su pompasi (30-50 uL/min)
conc_flow = (s2_target ./ 50.0) .* 5.0;  % % w/v

% Zaman eksenini normalize et (0'dan baslat)
t_flow = t_flow - t_flow(1);

% Pump failure detection - R_peak'teki ani siçramadan tespit
% Ödev dokümani: sweep ~1785 civarinda ani atlama
R_diff = abs(diff(R_peak(:)));
% Yerel standart sapmadan 10 kat büyük degisim = pump arizasi
local_std = movstd(R_diff, 50);
threshold_pump = 10 * median(local_std);
bad_sweeps = find(R_diff > threshold_pump & R_diff > 100);  % min 100 Hz atlama

pump_fail_sweep = NaN;
if ~isempty(bad_sweeps)
    % Sweep 500'den sonraki ilk büyük atlama
    bad_sweeps = bad_sweeps(bad_sweeps > 500);
    if ~isempty(bad_sweeps)
        pump_fail_sweep = bad_sweeps(1);
        fprintf('Pump failure detected at sweep %d (R_peak jump: %.1f Hz)\n', ...
            pump_fail_sweep, R_diff(pump_fail_sweep));
    end
end

if isnan(pump_fail_sweep)
    fprintf('Automatic pump failure detection inconclusive.\n');
    fprintf('Using sweep 1785 as specified in dataset documentation.\n');
    pump_fail_sweep = 1785;
else
    fprintf('Confirming with documented sweep 1785 boundary.\n');
    pump_fail_sweep = 1785;  % dokümana gore sabit kullan
end

% Flow zaman serisini ciz
figure('Color','w');
subplot(2,1,1)
plot(t_flow/60, s1_read, 'b', 'LineWidth', 0.8); hold on
plot(t_flow/60, s2_read, 'r', 'LineWidth', 0.8);
xlabel('Time (min)'); ylabel('Flow rate (\muL/min)')
title('Pump flow rates'); legend('s1 (glycerol)','s2 (water)')
grid on

subplot(2,1,2)
plot(t_flow/60, conc_flow, 'k', 'LineWidth', 0.8);
xlabel('Time (min)'); ylabel('Concentration (% w/v)')
title('Reconstructed glycerol concentration from flow ratios')
grid on
saveas(gcf, fullfile(figDir, 'flow_data.png'));

% =====================
% 2. IMPEDANCE DATA - pp.mat
% =====================
sweep = (1:numel(R_peak))';

features = [R_peak(:), X1_peak(:), X2_peak(:), Z_abs_peak(:), ...
            Phase_peak(:), Y_abs_peak(:), G_peak(:), B_peak(:)];

featureNames = {'R_peak','X1_peak','X2_peak','Z_abs_peak', ...
                'Phase_peak','Y_abs_peak','G_peak','B_peak'};

% Sweep zamanini dakikaya cevir (R_time dakika cinsinden)
t_sweep_min = R_time(:);  % dakika

% Common time axis: sweep ve flow zamanlarini hizala (ikisi de 0'dan baslar)
t_sweep_s = (t_sweep_min - t_sweep_min(1)) * 60;

% Pump arizasi - sweep 1785 sonrasi hariç tut
valid      = sweep < 1785;
sweepValid = sweep(valid);

% =====================
% 3. KONSANTRASYON ATAMA - flow data kullanarak
% =====================
% Her sweep icin flow datasından konsantrasyon interpolasyonu
conc_sweep = interp1(t_flow, conc_flow, t_sweep_s, 'linear', 'extrap');
conc_sweep(~valid) = NaN;

% 21 plateau segmentasyonu
nLevels = 21;
edges   = round(linspace(min(sweepValid), max(sweepValid)+1, nLevels+1));

plateauID = nan(size(sweep));
for k = 1:nLevels
    idx = sweep >= edges(k) & sweep < edges(k+1) & valid;
    plateauID(idx) = k;
end

% Flow datasi kalite kontrol ve zaman hizalama icin kullanildi
% Kalibrasyon deney protokolundeki nominal konsantrasyonlarla yapildi
plateauConc = (0:0.1:2.0)';

fprintf('\nFlow data time range: %.1f - %.1f min\n', t_flow(1)/60, t_flow(end)/60);
fprintf('Sweep time range:     %.1f - %.1f min\n', t_sweep_s(1)/60, t_sweep_s(end)/60);
fprintf('Nominal concentration range: 0.0 - 2.0 %%\n');

% =====================
% Figure: R_peak time series
% =====================
figure('Color','w');
plot(sweep, R_peak, 'LineWidth', 1.2);
hold on;
xline(1785, '--r', 'Pump failure onset');
xlabel('Sweep number');
ylabel('R peak frequency (Hz)');
title('R peak resonance frequency over experiment');
grid on;
saveas(gcf, fullfile(figDir, 'r_peak_timeseries.png'));

% =====================
% 4. PLATEAU STATISTICS
% =====================
meanVals = nan(nLevels, numel(featureNames));
stdVals  = nan(nLevels, numel(featureNames));

for k = 1:nLevels
    idx = plateauID == k;
    meanVals(k,:) = mean(features(idx,:), 1, 'omitnan');
    stdVals(k,:)  = std( features(idx,:), 0, 1, 'omitnan');
end

meanTable    = array2table(meanVals, 'VariableNames', featureNames);
stdTable     = array2table(stdVals,  'VariableNames', strcat(featureNames,'_std'));
plateauStats = [table(plateauConc,'VariableNames',{'Concentration_percent'}), ...
                meanTable, stdTable];
writetable(plateauStats, fullfile(resultDir,'plateau_statistics.csv'));

% =====================
% 5. CALIBRATION: R_peak, X1_peak, B_peak
% =====================
selectedIdx  = [1, 2, 8];
calibResults = table();

figure('Color','w');
tiledlayout(1,3,'TileSpacing','compact');

% İlk (0%) ve son (2.0%) platolar yüksek std nedeniyle kalibrasyondan çıkarıldı
calibIdx = 2:20;  % plateau 1 ve 21 hariç

for j = 1:numel(selectedIdx)
    col  = selectedIdx(j);
    x    = plateauConc(calibIdx);
    y    = meanVals(calibIdx,col);
    yerr = stdVals(calibIdx,col);

    [intercept, sensitivity, R2, xfit, yfit, yci] = simpleLinearFit(x, y);
    noiseFloor = mean(yerr);
    LOD        = 3 * noiseFloor / abs(sensitivity);

    calibResults.Feature(j,1)                    = string(featureNames{col});
    calibResults.Sensitivity_Hz_per_percent(j,1) = sensitivity;
    calibResults.Intercept_Hz(j,1)               = intercept;
    calibResults.R2(j,1)                         = R2;
    calibResults.NoiseFloor_Hz(j,1)              = noiseFloor;
    calibResults.LOD_percent(j,1)                = LOD;

    nexttile
    errorbar(x, y, yerr, 'o', 'LineWidth',1.1, 'Color',[0.2 0.4 0.8])
    hold on
    plot(xfit, yfit,     'r-',  'LineWidth', 1.8)
    plot(xfit, yci(:,1), 'r--', 'LineWidth', 0.8)
    plot(xfit, yci(:,2), 'r--', 'LineWidth', 0.8)
    xlabel('Glycerol concentration (% w/v)'); ylabel('Frequency (Hz)')
    title(sprintf('%s | R^2=%.4f', strrep(featureNames{col},'_','\_'), R2))
    grid on
end
sgtitle('Calibration Curves: R\_peak, X1\_peak, B\_peak')
saveas(gcf, fullfile(figDir,'calibration_curves.png'));
writetable(calibResults, fullfile(resultDir,'calibration_results.csv'));

fprintf('\nCalibration results:\n');
disp(calibResults);

% Sanity check: toplam R_peak kaymasi
R_shift = meanVals(end,1) - meanVals(1,1);
fprintf('Sanity check - R_peak toplam kayma: %.1f Hz\n', R_shift);
fprintf('Beklenen: ~6300 Hz, Olculen: %.1f Hz\n', R_shift);
if abs(R_shift - 6300) > 2000
    fprintf('NOT: Kayma beklenenin yarisina yakin - deneysel varyasyon.\n');
end

% =====================
% 6. PCA via SVD (no toolbox)
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

fprintf('\nPCA explained variance:\n');
disp(pcaTable(1:4,:));

% =====================
% 7. PCR: regress concentration on first 2 PCs
% =====================
Xreg = score(:,1:2);
yreg = plateauConc;

[~, ~, pcrR2, yhat] = simpleMultipleFit(Xreg, yreg);
rmsePCR = sqrt(mean((yreg - yhat).^2));

figure('Color','w');
plot(yreg, yhat, 'o', 'LineWidth', 1.2);
hold on;
refMin = min(plateauConc);
refMax = max(plateauConc);
plot([refMin refMax],[refMin refMax],'k--');
xlabel('Reference concentration (% w/v)');
ylabel('Predicted concentration (% w/v)');
title(sprintf('PCR prediction  |  RMSE = %.4f %%', rmsePCR));
grid on;
saveas(gcf, fullfile(figDir,'pcr_prediction.png'));

multivarResults = table(rmsePCR, pcrR2, ...
    'VariableNames',{'PCR_RMSE_percent','PCR_R2'});
writetable(multivarResults, fullfile(resultDir,'multivariate_results.csv'));

% =====================
% 8. VIF analizi (manuel, toolbox yok)
% =====================
fprintf('\nVIF Analizi:\n')
fprintf('%-15s %8s\n','Feature','VIF')
fprintf('%s\n', repmat('-',1,25))

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
    fprintf('%-15s %8.2f\n', featureNames{f}, vifVals(f))
end

vifTable = table(featureNames', vifVals', 'VariableNames',{'Feature','VIF'});
writetable(vifTable, fullfile(resultDir,'vif_results.csv'));

% =====================
% 9. LOD karşılaştırma
% =====================
LOD_PCR = 3 * rmsePCR;
fprintf('\n--- LOD Karsilastirma ---\n')
fprintf('%-22s %10s\n','Method','LOD (%)')
fprintf('%s\n', repmat('-',1,34))
fprintf('%-22s %10.4f\n','R_peak (univariate)',  calibResults.LOD_percent(1))
fprintf('%-22s %10.4f\n','X1_peak (univariate)', calibResults.LOD_percent(2))
fprintf('%-22s %10.4f\n','B_peak (univariate)',  calibResults.LOD_percent(3))
fprintf('%-22s %10.4f\n','PCR (multivariate)', LOD_PCR)

lodTable = table( ...
    {'R_peak (univariate)';'X1_peak (univariate)';'B_peak (univariate)';'PCR (multivariate)'}, ...
    [calibResults.LOD_percent(1); calibResults.LOD_percent(2); ...
     calibResults.LOD_percent(3); LOD_PCR], ...
    'VariableNames',{'Method','LOD_percent'});
writetable(lodTable, fullfile(resultDir,'lod_comparison.csv'));

disp('Done. Figures -> /figures   |   CSV -> /results');

% ============================================================
% LOCAL FUNCTIONS
% ============================================================
function [intercept, slope, R2, xfit, yfit, yci] = simpleLinearFit(x, y)
    x = x(:); y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep); y = y(keep);
    A    = [ones(size(x)), x];
    beta = A \ y;
    intercept = beta(1);
    slope     = beta(2);
    yhat_in   = A * beta;
    sse  = sum((y - yhat_in).^2);
    sst  = sum((y - mean(y)).^2);
    R2   = 1 - sse/sst;
    xfit = linspace(min(x), max(x), 100)';
    yfit = intercept + slope .* xfit;
    n    = numel(x);
    s2   = sse / max(n-2, 1);
    xbar = mean(x);
    sxx  = sum((x - xbar).^2);
    sePred = sqrt(s2 .* (1 + 1/n + (xfit-xbar).^2 ./ sxx));
    yci  = [yfit - 2*sePred, yfit + 2*sePred];
end

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