%% =========================================================
%  mdpi_poly_data4.m
%  data4 ham spektrum dosyalarindan 3. dereceden polinom fit
%  Cikti: 3rd_fit_data4.mat (T tablosu)
%
%  Calistirmadan once MATLAB'i data4 klasorune almayi unutma:
%  cd('....\data\data4')
% =========================================================
clearvars
clc

%% --- Dosya yolunu ayarla ---
data_path = '.\';   % data4 klasoru

%% --- Veriyi yukle ---
% Her dosya: [nSweep x 3000] boyutunda
% Yeniden sekillendirme: [nSweep x 1000 x 3]
%   3. boyut: (1) frekans ekseni, (2) real kisim, (3) imag kisim

fprintf('Dosyalar yukleniyor...\n');

raw = readmatrix(fullfile(data_path, 'B1_data.txt'));
B1_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'B2_data.txt'));
B2_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'angle_Z_data.txt'));
angle_Z_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'Z_abs_data.txt'));
Z_abs_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'Z_abs2_data.txt'));
Z_abs2_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'X1_data.txt'));
X1_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'X2_data.txt'));
X2_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'R_data.txt'));
R_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

raw = readmatrix(fullfile(data_path, 'G_data.txt'));
G_data = permute(reshape(raw.', 3, 1000, []), [3 2 1]);

nSweep = size(R_data, 1);
fprintf('  Toplam sweep: %d\n', nSweep);

%% --- Her sweep icin polinom fit ---
fprintf('Polinom fit yapiliyor...\n');

B1_poly      = zeros(nSweep, 4);
B2_poly      = zeros(nSweep, 4);
angle_Z_poly = zeros(nSweep, 4);
Z_abs_poly   = zeros(nSweep, 4);
Z_abs2_poly  = zeros(nSweep, 4);
X1_poly      = zeros(nSweep, 4);
X2_poly      = zeros(nSweep, 4);
R_poly       = zeros(nSweep, 4);

for k = 1:nSweep
    B1_poly(k,:)      = mdpi_3rd_poly([B1_data(k,:,1)',      imag(1./(B1_data(k,:,2)      + 1i*B1_data(k,:,3)))']);
    B2_poly(k,:)      = mdpi_3rd_poly([B2_data(k,:,1)',      imag(1./(B2_data(k,:,2)      + 1i*B2_data(k,:,3)))']);
    angle_Z_poly(k,:) = mdpi_3rd_poly([angle_Z_data(k,:,1)', atan(angle_Z_data(k,:,3) ./ angle_Z_data(k,:,2))']);
    Z_abs_poly(k,:)   = mdpi_3rd_poly([Z_abs_data(k,:,1)',   sqrt(Z_abs_data(k,:,2).^2  + Z_abs_data(k,:,3).^2)']);
    Z_abs2_poly(k,:)  = mdpi_3rd_poly([Z_abs2_data(k,:,1)',  sqrt(Z_abs2_data(k,:,2).^2 + Z_abs2_data(k,:,3).^2)']);
    X1_poly(k,:)      = mdpi_3rd_poly([X1_data(k,:,1)',      X1_data(k,:,3)']);
    X2_poly(k,:)      = mdpi_3rd_poly([X2_data(k,:,1)',      X2_data(k,:,3)']);
    R_poly(k,:)       = mdpi_3rd_poly([R_data(k,:,1)',       R_data(k,:,2)']);
end

fprintf('  Fit tamamlandi.\n');

%% --- Tabloyu olustur ve kaydet ---
T = array2table([B1_poly(:,1:4), B2_poly(:,1:4), angle_Z_poly(:,1:4), ...
                 Z_abs_poly(:,1:4), Z_abs2_poly(:,1:4), ...
                 X1_poly(:,1:4), X2_poly(:,1:4), R_poly(:,1:4)], ...
    'VariableNames', { ...
        'B1_a','B1_b','B1_c','B1_d', ...
        'B2_a','B2_b','B2_c','B2_d', ...
        'angle_Z_a','angle_Z_b','angle_Z_c','angle_Z_d', ...
        'Z_abs_a','Z_abs_b','Z_abs_c','Z_abs_d', ...
        'Z_abs2_a','Z_abs2_b','Z_abs2_c','Z_abs2_d', ...
        'X1_a','X1_b','X1_c','X1_d', ...
        'X2_a','X2_b','X2_c','X2_d', ...
        'R_a','R_b','R_c','R_d'});

save('3rd_fit_data4.mat', 'T');
fprintf('  3rd_fit_data4.mat kaydedildi.\n');

head(T)
summary(T)