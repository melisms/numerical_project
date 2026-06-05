clearvars
clc
dataDir = fileparts(mfilename('fullpath'));
dataRoot = fileparts(dataDir);
%% =========================
% 1-2. VERİ & ÖZELLİK MATRİSİ
% =========================
if isfile(fullfile(dataRoot,'3rd_fit_v3.mat'))
    fprintf('3rd_fit_v3.mat yükleniyor...\n');
    load('3rd_fit_v3.mat','T_all');
    X_all     = table2array(T_all);
    all_names = T_all.Properties.VariableNames;
else
    error('Önce feature_engineering_v3.m çalıştırın.');
end

%% =========================
% 3. HEDEF DEĞİŞKEN
% =========================
pp = load("pp.mat");
baseDir2 = fileparts(dataDir);
dataFile = fullfile(baseDir2, 'IWIS2025flow', 'data.txt');
opts = detectImportOptions(dataFile);
opts.VariableNamingRule = 'modify';
Traw = readtable(dataFile, opts);
s1=Traw{:,4}; s2=Traw{:,5};
conc=(s1./(s1+s2))*2.0;
t_sweep_s=pp.R_time(:)*60;
t_flow=Traw{:,2}; t_flow=t_flow-t_flow(1);
y=interp1(t_flow, conc, t_sweep_s, 'linear', 'extrap'); y=y(:);

%% =========================
% 4. TEMİZLİK & STANDARDİZASYON
% =========================
bad_cols = any(isnan(X_all)|isinf(X_all),1) | (std(X_all)==0);
X_all(:,bad_cols)=[]; all_names(bad_cols)=[];
bad_rows = any(isnan(X_all)|isinf(X_all),2)|isnan(y);
X_all(bad_rows,:)=[]; y(bad_rows)=[];
fprintf('Temizlenen: %d sütun, %d satır\n', sum(bad_cols), sum(bad_rows));

mu_x=mean(X_all); sig_x=std(X_all); sig_x(sig_x==0)=1;
Xz=(X_all-mu_x)./sig_x;

%% =========================
% 5. TRAIN / TEST SPLIT
% =========================
rng(42);
n=size(Xz,1);
idx=randperm(n);
nTrain=round(0.7*n);
trainIdx=idx(1:nTrain); testIdx=idx(nTrain+1:end);
Xtrain=Xz(trainIdx,:); ytrain=y(trainIdx);
Xtest=Xz(testIdx,:);   ytest=y(testIdx);

R2fn   = @(yt,yp) 1-sum((yt-yp).^2)/sum((yt-mean(yt)).^2);
RMSEfn = @(yt,yp) sqrt(mean((yt-yp).^2));

%% =========================
% 6. ÖZELLİK SEÇİMİ — Permutation Importance
% =========================
fprintf('\nPermutation importance hesaplanıyor...\n');
rng(42);
mdl_base = fitrensemble(Xtrain, ytrain,'Method','Bag', ...
    'NumLearningCycles',200,'Learners',templateTree('MinLeafSize',3));

yp_base  = predict(mdl_base, Xtest);
r2_base  = R2fn(ytest, yp_base);
fprintf('Baseline R² (tüm özellikler): %.4f\n', r2_base);

nFeat = size(Xtest,2);
perm_imp = zeros(1,nFeat);
rng(0);
for fi = 1:nFeat
    Xtest_perm = Xtest;
    Xtest_perm(:,fi) = Xtest_perm(randperm(size(Xtest,1)),fi);
    yp_perm = predict(mdl_base, Xtest_perm);
    perm_imp(fi) = r2_base - R2fn(ytest, yp_perm);  
end

[~, perm_ord] = sort(perm_imp,'descend');

fprintf('\nEn önemli 15 özellik (permutation):\n');
for i = 1:15
    fprintf('  %2d. %-26s Δ R²=%.4f\n', i, all_names{perm_ord(i)}, perm_imp(perm_ord(i)));
end

%% =========================
% 7. TOP-N TARAMA — Permutation sıralamasıyla
% =========================
top_ns = [5,10,15,20,30,50,75];
r2_scan = zeros(size(top_ns));
rmse_scan = zeros(size(top_ns));

fprintf('\n%-8s  %8s  %8s\n','Top-N','R²','RMSE');
fprintf('----------------------------\n');
for ni = 1:numel(top_ns)
    topN = min(top_ns(ni), nFeat);
    sel  = perm_ord(1:topN);
    rng(42);
    m = fitrensemble(Xtrain(:,sel), ytrain,'Method','Bag', ...
        'NumLearningCycles',300,'Learners',templateTree('MinLeafSize',3));
    yp = predict(m, Xtest(:,sel));
    r2_scan(ni)   = R2fn(ytest,yp);
    rmse_scan(ni) = RMSEfn(ytest,yp);
    fprintf('%-8d  %8.4f  %8.4f\n', topN, r2_scan(ni), rmse_scan(ni));
end

[~,best_ni] = max(r2_scan);
best_N = top_ns(best_ni);
sel_final = perm_ord(1:best_N);
fprintf('\nEn iyi Top-N: %d  (R²=%.4f)\n', best_N, r2_scan(best_ni));

%% =========================
% 8. FİNAL MODEL
% =========================
fprintf('\nFinal model eğitiliyor...\n');
rng(42);
mdl_final = fitrensemble(Xtrain(:,sel_final), ytrain,'Method','Bag', ...
    'NumLearningCycles',500,'Learners', ...
    templateTree('MinLeafSize',3,'MaxNumSplits',60));

% Train & Test performansı
yp_train = predict(mdl_final, Xtrain(:,sel_final));
y_pred   = predict(mdl_final, Xtest(:,sel_final));

fprintf('\n=========================================\n');
fprintf('FİNAL MODEL PERFORMANSI\n');
fprintf('-----------------------------------------\n');
fprintf('Train R²  = %.4f  (overfit kontrolü)\n', R2fn(ytrain, yp_train));
fprintf('Test  R²  = %.4f\n', R2fn(ytest, y_pred));
fprintf('Test RMSE = %.4f\n', RMSEfn(ytest, y_pred));
fprintf('Test MAE  = %.4f\n', mean(abs(ytest-y_pred)));
fprintf('Özellik sayısı: %d\n', best_N);
fprintf('=========================================\n');

%% =========================
% 9. GRAFİKLER
% =========================
figure('Name','Final Model','Position',[50 400 1100 380]);
plot(ytest,'b-','LineWidth',1.5,'DisplayName','Actual'); hold on;
plot(y_pred,'r--','LineWidth',1.5,'DisplayName','Predicted');
legend; grid on;
title(sprintf('Final RF (Test R²=%.4f, RMSE=%.4f)', R2fn(ytest,y_pred), RMSEfn(ytest,y_pred)));
xlabel('Test Sample'); ylabel('Concentration');

figure('Name','Parity Plot','Position',[50 50 460 420]);
mn=min([ytest;y_pred]); mx=max([ytest;y_pred]);
scatter(ytest,y_pred,40,'filled','MarkerFaceAlpha',0.6); hold on;
plot([mn mx],[mn mx],'r-','LineWidth',1.5);
xlabel('Actual'); ylabel('Predicted'); grid on; axis equal;
title(sprintf('Parity Plot (R²=%.4f)', R2fn(ytest,y_pred)));

figure('Name','Residual','Position',[560 50 680 380]);
res=ytest-y_pred;
subplot(1,2,1); plot(res,'Color',[0.2 0.5 0.8]); yline(0,'r--');
title('Residuals'); xlabel('Sample'); ylabel('Error'); grid on;
subplot(1,2,2); histogram(res,25,'FaceColor',[0.2 0.5 0.8],'EdgeColor','w');
title('Distribution'); xlabel('Error'); ylabel('Frequency'); grid on;

figure('Name','Permutation Importance','Position',[50 50 800 380]);
top_show = min(20, nFeat);
barh(flip(perm_imp(perm_ord(1:top_show))),'FaceColor',[0.3 0.65 0.4]);
yticks(1:top_show);
yticklabels(flip(all_names(perm_ord(1:top_show))));
title('Top-20 Features (Permutation Importance)'); xlabel('ΔR² (loss when feature is removed)'); grid on;

figure('Name','Top-N Scan','Position',[50 50 500 340]);
plot(top_ns, r2_scan,'b-o','LineWidth',1.5,'MarkerSize',7);
xline(best_N,'r--',sprintf('Best: %d',best_N));
xlabel('Number of Top-N Features'); ylabel('Test R²'); title('Feature Count Scan'); grid on;