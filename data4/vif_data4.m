%% =========================================================
%  vif_data4.m
%  data4 VIF (Varyans Enflasyon Faktoru) Analizi
%  Ozellikler arasi coklu dogrusalligi olcer
%
%  On kosul: 3rd_fit_data4.mat mevcut olmali
%            (mdpi_poly_data4.m calistirilmis olmali)
%  KULLANIM: data4 klasorunde calistir
% =========================================================
clearvars
clc
close all

%% =========================
% 1. VERI YUKLE
% =========================
load('3rd_fit_data4.mat', 'T');
X = table2array(T);
featureNames = T.Properties.VariableNames;

fprintf('Yuklendi: %d sweep x %d ozellik\n', size(X,1), size(X,2));

%% =========================
% 2. SABiT OZELLIK CIKAR
% =========================
sigma0 = std(X);
const_idx = (sigma0 == 0);
if any(const_idx)
    fprintf('Sabit ozellik cikariliyor: %d adet\n', sum(const_idx));
    X(:, const_idx) = [];
    featureNames(const_idx) = [];
end
nFeat = size(X, 2);

%% =========================
% 3. STANDARDiZASYON
% =========================
mu    = mean(X, 1);
sigma = std(X, 0, 1);
sigma(sigma == 0) = 1;
Xz = (X - mu) ./ sigma;

%% =========================
% 4. VIF HESABI
%    Her ozellik icin: digerlerinden o ozelligi tahmin et,
%    VIF = 1 / (1 - R2)
% =========================
fprintf('\n=== VIF ANALiZi (data4) ===\n');
fprintf('%-14s %14s\n', 'Feature', 'VIF');
fprintf('%s\n', repmat('-', 1, 30));

vifVals = zeros(1, nFeat);
for f = 1:nFeat
    y_vif  = Xz(:, f);
    X_vif  = Xz(:, setdiff(1:nFeat, f));

    % OLS ile R2
    b_vif  = X_vif \ y_vif;
    yh_vif = X_vif * b_vif;
    ss_res = sum((y_vif - yh_vif).^2);
    ss_tot = sum((y_vif - mean(y_vif)).^2);
    R2_vif = 1 - ss_res / ss_tot;

    vifVals(f) = 1 / max(1 - R2_vif, 1e-10);
    fprintf('%-14s %14.2f\n', featureNames{f}, vifVals(f));
end

%% =========================
% 5. OZET ISTATISTIKLER
% =========================
high_vif = sum(vifVals > 10);
[maxVIF, maxIdx] = max(vifVals);   % sadece ilk maksimumu al
fprintf('%s\n', repmat('-', 1, 30));
fprintf('VIF > 10 olan ozellik sayisi: %d / %d\n', high_vif, nFeat);
fprintf('Medyan VIF: %.2e\n', median(vifVals));
fprintf('Maksimum VIF: %.2e (%s)\n', maxVIF, featureNames{maxIdx});
fprintf('NOT: VIF>10^5 degerleri asiri collinearity gosterir.\n');

%% =========================
% 6. CSV KAYDET
% =========================
if ~exist('results', 'dir'), mkdir('results'); end
vifTable = table(featureNames', vifVals', ...
    'VariableNames', {'Feature', 'VIF'});
writetable(vifTable, fullfile('results', 'vif_results_data4.csv'));
fprintf('\nresults/vif_results_data4.csv kaydedildi.\n');

%% =========================
% 7. GRAFIK
% =========================
figure('Name', 'VIF Analizi', 'Position', [100 100 900 450]);
bar(vifVals, 'FaceColor', [0.4 0.7 0.4]);
hold on;
yline(10, 'r--', 'Esik = 10', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');   % VIF cok genis araliktaysa log olcek
xticks(1:nFeat);
xticklabels(featureNames);
xtickangle(45);
ylabel('VIF Degeri (log olcek)');
title('Varyans Enflasyon Faktoru (VIF) — data4');
grid on;
saveas(gcf, fullfile('results', 'vif_plot_data4.png'));

fprintf('VIF grafigi kaydedildi.\n');