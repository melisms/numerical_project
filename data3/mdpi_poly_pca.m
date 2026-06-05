clearvars
clc

dataDir = fileparts(mfilename('fullpath'));
% Load the datasets
B1_data      = loadDataset(fullfile(dataDir,'B1_data.txt'));
B2_data      = loadDataset(fullfile(dataDir,'B2_data.txt'));
angle_Z_data = loadDataset(fullfile(dataDir,'angle_Z_data.txt'));
Z_abs_data   = loadDataset(fullfile(dataDir,'Z_abs_data.txt'));
Z_abs2_data  = loadDataset(fullfile(dataDir,'Z_abs2_data.txt'));
X1_data      = loadDataset(fullfile(dataDir,'X1_data.txt'));
X2_data      = loadDataset(fullfile(dataDir,'X2_data.txt'));
R_data       = loadDataset(fullfile(dataDir,'R_data.txt'));
G_data       = loadDataset(fullfile(dataDir,'G_data.txt'));

for k=1:size(B1_data,1)
    B1_poly(k,:)=mdpi_3rd_poly([B1_data(k,:,1)',imag(1./(B1_data(k,:,2)+i*B1_data(k,:,3)))']);
    B2_poly(k,:)=mdpi_3rd_poly([B2_data(k,:,1)',imag(1./(B2_data(k,:,2)+i*B2_data(k,:,3)))']);
    angle_Z_poly(k,:)=mdpi_3rd_poly([angle_Z_data(k,:,1)',atan(angle_Z_data(k,:,3)./angle_Z_data(k,:,2))']);
    Z_abs_poly(k,:)=mdpi_3rd_poly([Z_abs_data(k,:,1)',sqrt(Z_abs_data(k,:,2).^2+Z_abs_data(k,:,3).^2)']);
    Z_abs2_poly(k,:)=mdpi_3rd_poly([Z_abs2_data(k,:,1)',sqrt(Z_abs2_data(k,:,2).^2+Z_abs2_data(k,:,3).^2)']);
    X1_poly(k,:)=mdpi_3rd_poly([X1_data(k,:,1)',X1_data(k,:,3)']);
    X2_poly(k,:)=mdpi_3rd_poly([X2_data(k,:,1)',X2_data(k,:,3)']);
    R_poly(k,:)=mdpi_3rd_poly([R_data(k,:,1)',R_data(k,:,2)']);
    % G_poly(k,:)=mdpi_3rd_poly([G_data(k,:,1)',real(1./(G_data(k,:,2)+i*G_data(k,:,3)))']);
end
T = array2table([B1_poly(:,1:4),B2_poly(:,1:4),angle_Z_poly(:,1:4),Z_abs_poly(:,1:4),Z_abs2_poly(:,1:4),X1_poly(:,1:4),X2_poly(:,1:4),R_poly(:,1:4)],...
    'VariableNames',{'B1_a','B1_b','B1_c','B1_d','B2_a','B2_b','B2_c','B2_d','angle_Z_a','angle_Z_b','angle_Z_c','angle_Z_d',...
    'Z_abs_a','Z_abs_b','Z_abs_c','Z_abs_d','Z_abs2_a','Z_abs2_b','Z_abs2_c','Z_abs2_d','X1_a','X1_b','X1_c','X1_d',...
    'X2_a','X2_b','X2_c','X2_d','R_a','R_b','R_c','R_d'});
save('3rd_fit.mat',"T");

head(T)
summary(T)
X = table2array(T);

function out = loadDataset(filename)

    raw = readmatrix(filename);
    
    out = permute( ...
        reshape(raw.',3,1000,[]), ...
        [3 2 1]);

end