function [distortion] = THD(ak,bk,N)

c_1 = (((ak(1))^2 + (bk(1))^2)^(1/2)) / (2^(1/2));

sum_n = 0;

for k = 2:N
    c_n(k) = ((((ak(k))^2 + (bk(k))^2)^(1/2)) / (2^(1/2)))^2;
    sum_n = c_n(k) + sum_n;
end

sum_n = sqrt(sum_n);

distortion = sum_n / c_1;