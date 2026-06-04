%% =========================================================
%  stability_data4.m
%  data4 Kararlilik / Drift Analizi
%  On kosul: 3rd_fit_data4.mat ve pp.mat mevcut olmali
% =========================================================
clearvars
clc
close all

%% =========================
% 1. VERI YUKLE
% =========================
load('3rd_fit_data4.mat', 'T');
pp = load('pp.mat');

time_min = pp.R_time(:);   % dakika
R_peak   = pp.R_peak(:);
X1_peak  = pp.X1_peak(:);
X2_peak  = pp.X2_peak(:);

nSweep = length(time_min);
fprintf('Toplam sweep: %d\n', nSweep);
fprintf('Sure: %.1f - %.1f dakika\n', min(time_min), max(time_min));

%% =========================
% 2. POMPA ARIZASINI TESPIT ET
%    (~87. dakikadaki ani spike)
% =========================
R_diff     = abs(diff(R_peak));
threshold  = mean(R_diff) + 5 * std(R_diff);
fault_idx  = find(R_diff > threshold, 1, 'first');

if ~isempty(fault_idx)
    fprintf('Pompa arizasi tespit edildi: sweep %d (%.1f dk)\n', ...
        fault_idx, time_min(fault_idx));
    stable_idx = 1:fault_idx-1;
else
    fprintf('Pompa arizasi tespit edilmedi.\n');
    stable_idx = 1:nSweep;
end

% Kararli bolge
t_stable  = time_min(stable_idx);
R_stable  = R_peak(stable_idx);
X1_stable = X1_peak(stable_idx);
X2_stable = X2_peak(stable_idx);

fprintf('Kararli bolge: %d sweep (%.1f dk)\n', ...
    length(stable_idx), max(t_stable));

%% =========================
% 3. GURULTU TABANI (NOISE FLOOR)
% =========================
% Kararli bolgedeki standart sapma = gurultu tabani
noise_R  = std(R_stable);
noise_X1 = std(X1_stable);
noise_X2 = std(X2_stable);

fprintf('\n=== GURULTU TABANI ===\n');
fprintf('  R_peak  std: %.2f Hz\n', noise_R);
fprintf('  X1_peak std: %.2f Hz\n', noise_X1);
fprintf('  X2_peak std: %.2f Hz\n', noise_X2);

%% =========================
% 4. DRIFT ANALiZi
%    Lineer fit ile uzun vadeli kaymayi hesapla
% =========================
p_R  = polyfit(t_stable, R_stable,  1);
p_X1 = polyfit(t_stable, X1_stable, 1);
p_X2 = polyfit(t_stable, X2_stable, 1);

% Drift Hz/dakika ve Hz/saat cinsinden
drift_R_per_min  = p_R(1);
drift_X1_per_min = p_X1(1);
drift_X2_per_min = p_X2(1);

fprintf('\n=== DRiFT ANALiZi ===\n');
fprintf('  R_peak  drift: %.4f Hz/dk  (%.2f Hz/saat)\n', ...
    drift_R_per_min,  drift_R_per_min*60);
fprintf('  X1_peak drift: %.4f Hz/dk  (%.2f Hz/saat)\n', ...
    drift_X1_per_min, drift_X1_per_min*60);
fprintf('  X2_peak drift: %.4f Hz/dk  (%.2f Hz/saat)\n', ...
    drift_X2_per_min, drift_X2_per_min*60);

%% =========================
% 5. LOD HESABI (kararlilik bazli)
%    LOD = 3 * noise / sensitivity
%    sensitivity: data3 referans degeri ~6300 Hz / 2% = 3150 Hz/%
% =========================
sensitivity_ref = 3150;   % Hz/% (data3'ten beklenen)

LOD_R  = 3 * noise_R  / sensitivity_ref;
LOD_X1 = 3 * noise_X1 / sensitivity_ref;
LOD_X2 = 3 * noise_X2 / sensitivity_ref;

fprintf('\n=== LOD (kararlilik bazli) ===\n');
fprintf('  R_peak  LOD: %.5f %%\n', LOD_R);
fprintf('  X1_peak LOD: %.5f %%\n', LOD_X1);
fprintf('  X2_peak LOD: %.5f %%\n', LOD_X2);

%% =========================
% 6. GRAFIK 1: Zaman Serisi + Pompa Arizasi
% =========================
figure('Name','Zaman Serisi','Position',[100 100 1000 450]);

subplot(2,1,1);
plot(time_min, R_peak - R_peak(1), 'b-', 'LineWidth', 1.0);
hold on;
if ~isempty(fault_idx)
    xline(time_min(fault_idx), 'r--', 'LineWidth', 2, ...
        'Label', 'Pompa Arizasi');
end
ylabel('\DeltaR_{peak} (Hz)');
title('data4 — R_{peak} Zaman Serisi');
grid on;

subplot(2,1,2);
plot(time_min, X1_peak - X1_peak(1), 'Color',[0.1 0.6 0.1], 'LineWidth', 1.0);
hold on;
if ~isempty(fault_idx)
    xline(time_min(fault_idx), 'r--', 'LineWidth', 2);
end
xlabel('Zaman (dakika)');
ylabel('\DeltaX1_{peak} (Hz)');
title('data4 — X1_{peak} Zaman Serisi');
grid on;

%% =========================
% 7. GRAFIK 2: Drift Analizi
% =========================
figure('Name','Drift Analizi','Position',[100 580 1000 400]);

subplot(1,2,1);
plot(t_stable, R_stable - R_stable(1), 'b.', 'MarkerSize', 3);
hold on;
fit_line = polyval(p_R, t_stable) - polyval(p_R, t_stable(1));
plot(t_stable, fit_line, 'r-', 'LineWidth', 2);
xlabel('Zaman (dakika)');
ylabel('\DeltaR_{peak} (Hz)');
title(sprintf('R_{peak} Drift\n%.4f Hz/dk (%.2f Hz/saat)', ...
    drift_R_per_min, drift_R_per_min*60));
legend('Veri', 'Lineer Fit');
grid on;

subplot(1,2,2);
plot(t_stable, X1_stable - X1_stable(1), 'g.', 'MarkerSize', 3);
hold on;
fit_line2 = polyval(p_X1, t_stable) - polyval(p_X1, t_stable(1));
plot(t_stable, fit_line2, 'r-', 'LineWidth', 2);
xlabel('Zaman (dakika)');
ylabel('\DeltaX1_{peak} (Hz)');
title(sprintf('X1_{peak} Drift\n%.4f Hz/dk (%.2f Hz/saat)', ...
    drift_X1_per_min, drift_X1_per_min*60));
legend('Veri', 'Lineer Fit');
grid on;

sgtitle('Drift Analizi — data4 (Kararli Bolge)');

%% =========================
% 8. GRAFIK 3: Gurultu Dagilimi (Histogram)
% =========================
figure('Name','Gurultu Dagilimi','Position',[100 100 1000 380]);

subplot(1,3,1);
R_detrended = R_stable - polyval(p_R, t_stable);
histogram(R_detrended, 30, 'FaceColor',[0.2 0.4 0.8], 'Normalization','pdf');
xlabel('\DeltaR_{peak} (Hz)');
ylabel('PDF');
title(sprintf('R_{peak} Gurultu\n\\sigma = %.2f Hz', noise_R));
grid on;

subplot(1,3,2);
X1_detrended = X1_stable - polyval(p_X1, t_stable);
histogram(X1_detrended, 30, 'FaceColor',[0.1 0.6 0.1], 'Normalization','pdf');
xlabel('\DeltaX1_{peak} (Hz)');
ylabel('PDF');
title(sprintf('X1_{peak} Gurultu\n\\sigma = %.2f Hz', noise_X1));
grid on;

subplot(1,3,3);
X2_detrended = X2_stable - polyval(p_X2, t_stable);
histogram(X2_detrended, 30, 'FaceColor',[0.8 0.4 0.1], 'Normalization','pdf');
xlabel('\DeltaX2_{peak} (Hz)');
ylabel('PDF');
title(sprintf('X2_{peak} Gurultu\n\\sigma = %.2f Hz', noise_X2));
grid on;

sgtitle('Gurultu Dagilimi — data4');

%% =========================
% 9. OZET TABLO
% =========================
fprintf('\n========================================\n');
fprintf('       KARARLILIK ANALiZi OZETI\n');
fprintf('========================================\n');
fprintf('  Kararli sure    : %.1f dakika\n', max(t_stable));
fprintf('  Toplam sweep    : %d\n', length(stable_idx));
fprintf('\n  Ozellik   | Noise(Hz) | Drift(Hz/saat) | LOD(%%)\n');
fprintf('  ----------|-----------|----------------|--------\n');
fprintf('  R_peak    | %9.2f | %14.2f | %.5f\n', ...
    noise_R,  drift_R_per_min*60,  LOD_R);
fprintf('  X1_peak   | %9.2f | %14.2f | %.5f\n', ...
    noise_X1, drift_X1_per_min*60, LOD_X1);
fprintf('  X2_peak   | %9.2f | %14.2f | %.5f\n', ...
    noise_X2, drift_X2_per_min*60, LOD_X2);
fprintf('========================================\n');