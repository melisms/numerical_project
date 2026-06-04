function p = mdpi_3rd_poly(data)
% mdpi_3rd_poly — Her sweep icin 3. dereceden polinom fit yapar
% Girdi: data = [frekans_vektoru, sinyal_vektoru] (Nx2 matris)
% Cikti: p = polinom katsayilari (1x4)
%
% data2, data3, data4 icin ortaktir — degistirme.

x = data(:,1);
y = data(:,2);
[p,~,~] = polyfit(x, y, 3);
end