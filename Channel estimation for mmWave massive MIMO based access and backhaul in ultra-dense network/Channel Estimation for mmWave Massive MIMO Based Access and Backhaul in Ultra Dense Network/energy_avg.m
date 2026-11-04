function z = energy_avg(M)
m = size(M,2); 
z = sqrt(sum(abs(M).^2,2)/m);