clearvars
clc
dataDir = fileparts(mfilename('fullpath'));

%% =========================
% 1. VERİ YÜKLEME
% =========================
B1_data      = loadDataset(fullfile(dataDir,'B1_data.txt'));
B2_data      = loadDataset(fullfile(dataDir,'B2_data.txt'));
angle_Z_data = loadDataset(fullfile(dataDir,'angle_Z_data.txt'));
Z_abs_data   = loadDataset(fullfile(dataDir,'Z_abs_data.txt'));
Z_abs2_data  = loadDataset(fullfile(dataDir,'Z_abs2_data.txt'));
X1_data      = loadDataset(fullfile(dataDir,'X1_data.txt'));
X2_data      = loadDataset(fullfile(dataDir,'X2_data.txt'));
R_data       = loadDataset(fullfile(dataDir,'R_data.txt'));
G_data       = loadDataset(fullfile(dataDir,'G_data.txt'));

nSweep = size(B1_data,1);

%% =========================
% 2. SİNYAL FONKSİYONLARI
% =========================
sig_fns  = { ...
    @(d) imag(1./(d(:,:,2)+1i*d(:,:,3))), ...   % B1
    @(d) imag(1./(d(:,:,2)+1i*d(:,:,3))), ...   % B2
    @(d) atan(d(:,:,3)./d(:,:,2)),         ...   % angle_Z
    @(d) sqrt(d(:,:,2).^2+d(:,:,3).^2),   ...   % Z_abs
    @(d) sqrt(d(:,:,2).^2+d(:,:,3).^2),   ...   % Z_abs2
    @(d) d(:,:,3),                         ...   % X1
    @(d) d(:,:,3),                         ...   % X2
    @(d) d(:,:,2),                         ...   % R
    @(d) real(1./(d(:,:,2)+1i*d(:,:,3)))   ...   % G
};
datasets = {B1_data,B2_data,angle_Z_data,Z_abs_data,Z_abs2_data,X1_data,X2_data,R_data,G_data};
ds_names = {'B1','B2','angle_Z','Z_abs','Z_abs2','X1','X2','R','G'};
nDS = numel(datasets);

%% =========================
% 3. POLİNOM KATSAYILARI
% =========================
poly_mat  = zeros(nSweep, nDS*4);
poly_names = {};
sfx = {'a','b','c','d'};
for k = 1:nSweep
    for di = 1:nDS
        d = datasets{di};
        freq_vec = d(k,:,1)';
        sig_vec  = sig_fns{di}(d(k,:,:))';
        coeffs   = mdpi_3rd_poly([freq_vec, sig_vec]);
        poly_mat(k,(di-1)*4+(1:4)) = coeffs(1:4);
    end
end
for di = 1:nDS
    for s = 1:4
        poly_names{end+1} = sprintf('%s_%s', ds_names{di}, sfx{s});
    end
end

%% =========================
% 4. İSTATİSTİKSEL ÖZELLİKLER
% =========================
stat_mat   = zeros(nSweep, nDS*7);
stat_names = {};
for di = 1:nDS
    sig = sig_fns{di}(datasets{di});
    col = (di-1)*7;
    stat_mat(:,col+1) = mean(sig,2);
    stat_mat(:,col+2) = std(sig,0,2);
    stat_mat(:,col+3) = max(sig,[],2);
    stat_mat(:,col+4) = min(sig,[],2);
    stat_mat(:,col+5) = max(sig,[],2)-min(sig,[],2);
    stat_mat(:,col+6) = mean(diff(sig,1,2),2);
    stat_mat(:,col+7) = trapz(sig,2);
    nm = ds_names{di};
    stat_names(end+(1:7)) = {[nm '_mean'],[nm '_std'],[nm '_max'], ...
                              [nm '_min'],[nm '_ptp'],[nm '_dfdx'],[nm '_area']};
end

%% =========================
% 5. YENİ: YOĞUN FREKANS ÖRNEKLEMESİ
%    - Alt band (1-500):   10 nokta (seyrek, düşük bilgi)
%    - Üst band (500-1000): 40 nokta (yoğun, yüksek bilgi)
% =========================
freq_idx_low  = round(linspace(1,   500,  10));
freq_idx_high = round(linspace(501, 1000, 40));
freq_idx_all  = unique([freq_idx_low, freq_idx_high]);

freq_mat   = zeros(nSweep, nDS*numel(freq_idx_all));
freq_names = {};
col = 0;
for di = 1:nDS
    sig = sig_fns{di}(datasets{di});
    for fi = 1:numel(freq_idx_all)
        col = col+1;
        freq_mat(:,col) = sig(:, freq_idx_all(fi));
        freq_names{end+1} = sprintf('%s_f%d', ds_names{di}, freq_idx_all(fi));
    end
end

%% =========================
% 6. YENİ: BANT İSTATİSTİKLERİ
%    Alt / Orta / Üst band için ayrı mean, std, area
% =========================
bands      = {1:333, 334:667, 668:1000};
band_names = {'low','mid','high'};
band_mat   = zeros(nSweep, nDS*numel(bands)*3);
bnd_names  = {};
col = 0;
for di = 1:nDS
    sig = sig_fns{di}(datasets{di});
    nm  = ds_names{di};
    for bi = 1:numel(bands)
        bsig = sig(:, bands{bi});
        col = col+1; band_mat(:,col) = mean(bsig,2);
        col = col+1; band_mat(:,col) = std(bsig,0,2);
        col = col+1; band_mat(:,col) = trapz(bsig,2);
        bn = band_names{bi};
        bnd_names(end+(1:3)) = {sprintf('%s_%s_mean',nm,bn), ...
                                  sprintf('%s_%s_std',nm,bn), ...
                                  sprintf('%s_%s_area',nm,bn)};
    end
end

%% =========================
% 7. YENİ: KANAL ORANLARI (ratio features)
%    Fiziksel olarak anlamlı oranlar
% =========================
Z_abs_sig  = sig_fns{4}(Z_abs_data);
Z_abs2_sig = sig_fns{5}(Z_abs2_data);
R_sig      = sig_fns{8}(R_data);
G_sig      = sig_fns{9}(G_data);
X1_sig     = sig_fns{6}(X1_data);

ratio_mat = zeros(nSweep, 5);
ratio_mat(:,1) = mean(Z_abs_sig,2)  ./ (mean(Z_abs2_sig,2)+1e-10);  % Z1/Z2
ratio_mat(:,2) = mean(R_sig,2)      ./ (mean(Z_abs_sig,2)+1e-10);   % R/|Z| = cos(phi)
ratio_mat(:,3) = mean(G_sig,2)      ./ (mean(R_sig,2)+1e-10);       % G/R
ratio_mat(:,4) = max(X1_sig,[],2)   ./ (mean(Z_abs_sig,2)+1e-10);   % X_peak/|Z|
ratio_mat(:,5) = std(Z_abs_sig,0,2) ./ (mean(Z_abs_sig,2)+1e-10);   % dispersiyon
ratio_names = {'Z1_Z2_ratio','cosine_phi','G_R_ratio','Xpeak_Z_ratio','Z_dispersion'};

%% =========================
% 8. TAM ÖZELLİK MATRİSİ
% =========================
X_all     = [poly_mat, stat_mat, freq_mat, band_mat, ratio_mat];
all_names = [poly_names, stat_names, freq_names, bnd_names, ratio_names];

T_all = array2table(X_all, 'VariableNames', all_names);
save('3rd_fit_v3.mat','T_all');

fprintf('=================================\n');
fprintf('Özellik Seti v3 Özeti\n');
fprintf('---------------------------------\n');
fprintf('Polinom katsayıları   : %d\n', size(poly_mat,2));
fprintf('İstatistiksel         : %d\n', size(stat_mat,2));
fprintf('Frekans (yoğun üst)   : %d\n', size(freq_mat,2));
fprintf('Bant istatistikleri   : %d\n', size(band_mat,2));
fprintf('Fiziksel oranlar      : %d\n', size(ratio_mat,2));
fprintf('---------------------------------\n');
fprintf('TOPLAM                : %d özellik\n', size(X_all,2));
fprintf('=================================\n');

%% =========================
% 9. KORELASYONu HESAPLA & KARŞILAŞTIR
% =========================
try
    pp = load("pp.mat");
    baseDir2 = fileparts(dataDir);
    dataFile = fullfile(baseDir2, 'IWIS2025flow', 'data.txt');
    opts = detectImportOptions(dataFile);
    opts = setvaropts(opts, opts.VariableNames{1}, 'InputFormat','MM/dd/uuuu hh:mm:ss aa');
    opts.VariableNamingRule = 'modify';
    Traw = readtable(dataFile, opts);
    s1 = Traw{:,4}; s2 = Traw{:,5};
    conc = (s1./(s1+s2))*2.0;
    t_sweep_s = pp.R_time(:)*60;
    t_flow = Traw{:,2}; t_flow = t_flow - t_flow(1);
    y = interp1(t_flow, conc, t_sweep_s, 'linear', 'extrap');
    y = y(:);

    Xz = (X_all - mean(X_all)) ./ (std(X_all)+1e-10);
    corr_all = abs(corr(Xz, y, 'rows','complete'));
    [sorted_c, order_c] = sort(corr_all,'descend');

    fprintf('\nEn yüksek 15 korelasyon (v3):\n');
    for i = 1:15
        fprintf('  %-26s r = %.4f\n', all_names{order_c(i)}, sorted_c(i));
    end
    fprintf('\nv2 max : 0.3566\n');
    fprintf('v3 max : %.4f\n', max(corr_all));
catch ME
    fprintf('Korelasyon hesaplanamadı: %s\n', ME.message);
end

%% =========================
% 10. MODEL (v3 özellik setiyle RF)
% =========================
try
    fprintf('\nRF modeli eğitiliyor (v3 özellik seti)...\n');
    rng(1);
    n = size(X_all,1);
    idx = randperm(n);
    nTrain = round(0.7*n);
    trainIdx = idx(1:nTrain); testIdx = idx(nTrain+1:end);

    mu_x = mean(X_all); sig_x = std(X_all); sig_x(sig_x==0) = 1;
    Xz_full = (X_all - mu_x) ./ sig_x;

    Xtrain = Xz_full(trainIdx,:); ytrain = y(trainIdx);
    Xtest  = Xz_full(testIdx,:);  ytest  = y(testIdx);

    % Özellik önemi ile seçim
    mdl_rf_full = fitrensemble(Xtrain, ytrain,'Method','Bag', ...
        'NumLearningCycles',300,'Learners',templateTree('MinLeafSize',3));
    imp = predictorImportance(mdl_rf_full);
    [~, imp_ord] = sort(imp,'descend');

    % Top-N tarama (10, 20, 30, 50)
    top_ns = [10, 20, 30, 50, 100];
    fprintf('\n%-10s  %8s  %8s\n','Top-N','R²','RMSE');
    fprintf('------------------------------\n');
    R2fn   = @(yt,yp) 1-sum((yt-yp).^2)/sum((yt-mean(yt)).^2);
    RMSEfn = @(yt,yp) sqrt(mean((yt-yp).^2));
    for ni = 1:numel(top_ns)
        topN = min(top_ns(ni), size(Xtrain,2));
        sel  = imp_ord(1:topN);
        m    = fitrensemble(Xtrain(:,sel), ytrain,'Method','Bag', ...
                   'NumLearningCycles',300,'Learners',templateTree('MinLeafSize',3));
        yp   = predict(m, Xtest(:,sel));
        fprintf('%-10d  %8.4f  %8.4f\n', topN, R2fn(ytest,yp), RMSEfn(ytest,yp));
    end
catch ME2
    fprintf('Model çalıştırılamadı: %s\n', ME2.message);
end

function out = loadDataset(filename)
    raw = readmatrix(filename);
    out = permute(reshape(raw.',3,1000,[]),[3 2 1]);
end