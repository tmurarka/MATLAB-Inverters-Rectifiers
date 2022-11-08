function [avg,ak,bk,rw,err] = fourser(t,x,T,N)
w_ac = 2*pi/T; % fundamental
delta_t = t(2)-t(1); % time step for fourier
period = round(T/delta_t); % data points going into one period

% loop for ak and bk
for k = 1:N
    ak(k) = 0;
    bk(k) = 0;
    for m = 1:period
        ak(k) = ak(k) + 2*(x(m)*cos(k*w_ac*t(m))*delta_t)/T; %finding ak value based on cos
        bk(k) = bk(k) + 2*(x(m)*sin(k*w_ac*t(m))*delta_t)/T; %finding bk value based on sin
    end
end

avg = avrg(x,T,delta_t); % average

rw = avg;


for k = 1:N
    rw = rw + ak(k)*cos(k*w_ac*t) + bk(k)*sin(k*w_ac*t); %fits fouirer coefficients based on cos
                                                         %and sin waves, additionally adds average
end

err = sqrt(abs(x-rw).^2); % error