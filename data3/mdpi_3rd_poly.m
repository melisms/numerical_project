function p = mdpi_3rd_poly(data)

x = data(:,1);
y = data(:,2);

[p,~,mu] = polyfit(x, y, 3);

p = p; % sadece coeff

end