function [state] = switching(D, time)

delta_t = 1/10000;
omega = 2*pi/delta_t;
stored_time = [0:delta_t/10:delta_t];
k = 1;
fourier = 0.5;

if time > delta_t
    interval = time - ((floor(time/delta_t))*delta_t);
else
    interval = time;
end

if interval < 0.00001
    index = 1;
elseif interval < 0.00002
    index = 2;
elseif interval < 0.00003
    index = 3;
elseif interval < 0.00004
    index = 4;
elseif interval < 0.00005
    index = 5;
elseif interval < 0.00006
    index = 6;
elseif interval < 0.00007
    index = 7;
elseif interval < 0.00008
    index = 8;
elseif interval < 0.00009
    index = 9;
else
    index = 10;
end

while(k < 10000)
    fourier = fourier + cos(k*omega*stored_time)*0.1*(-2/(k^2))*(1-cos(k*pi));
    k = k + 1;
end

if fourier(index) >= D
    state = 0;
else
    state = 1;
end

end
 