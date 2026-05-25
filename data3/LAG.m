clearvars
clc

%% =========================
% 1. LOAD DATA
% =========================
load('3rd_fit.mat','T');
X = table2array(T);

opts = detectImportOptions('/Users/melisamuslu/Documents/MATLAB/data/IWIS2025flow/data.txt');
opts.VariableNamingRule = 'modify';
Traw = readtable('/Users/melisamuslu/Documents/MATLAB/data/IWIS2025flow/data.txt', opts);

y = Traw{:,4};   % s1(Read)

%% =========================
% 2. ALIGN LENGTH
% =========================
n = min(size(X,1), length(y));
X = X(1:n,:);
y = y(1:n,:);

%% =========================
% 3. STANDARDIZATION
% =========================
mu = mean(X);
sigma = std(X);
sigma(sigma==0)=1;
Xz = (X - mu) ./ sigma;

%% =========================
% 4. FIND BEST LAG
% =========================
maxLag = 30;
bestR2 = -inf;
bestLag = 0;

for lag = 0:maxLag

    if lag == 0
        Xs = Xz;
        ys = y;
    else
        Xs = Xz(1:end-lag,:);
        ys = y(lag+1:end);
    end

    % train/test split (fixed per lag)
    idx = randperm(length(ys));
    nTrain = round(0.7*length(ys));

    tr = idx(1:nTrain);
    te = idx(nTrain+1:end);

    Xtr = Xs(tr,:);
    ytr = ys(tr);

    Xte = Xs(te,:);
    yte = ys(te);

    % SVR model
    mdl = fitrsvm(Xtr, ytr, ...
        'KernelFunction','gaussian', ...
        'Standardize',true);

    ypred = predict(mdl, Xte);

    % R2
    R2 = 1 - sum((yte - ypred).^2) / sum((yte - mean(yte)).^2);

    if R2 > bestR2
        bestR2 = R2;
        bestLag = lag;
        bestModel = mdl;
        bestTest = struct('Xte',Xte,'yte',yte,'ypred',ypred);
    end
end

fprintf('=========================\n');
fprintf('BEST LAG: %d\n', bestLag);
fprintf('BEST R2 : %.4f\n', bestR2);
fprintf('=========================\n');

%% =========================
% 5. FINAL EVALUATION
% =========================
ytest = bestTest.yte;
ypred = bestTest.ypred;

RMSE = sqrt(mean((ytest - ypred).^2));
R2 = bestR2;

fprintf('FINAL MODEL\n');
fprintf('R2   = %.4f\n', R2);
fprintf('RMSE = %.4f\n', RMSE);

%% =========================
% 6. PLOTS
% =========================
figure;
plot(ytest,'b','LineWidth',1.5); hold on;
plot(ypred,'r--','LineWidth',1.5);
legend('True','Predicted');
title(['Lag-Optimized SVR (Lag = ' num2str(bestLag) ')']);
grid on;

figure;
scatter(ytest, ypred, 25, 'filled');
xlabel('True');
ylabel('Predicted');
title('Parity Plot (Lag Optimized)');
grid on;
refline(1,0);

figure;
residuals = ytest - ypred;
plot(residuals);
title('Residuals');
grid on;