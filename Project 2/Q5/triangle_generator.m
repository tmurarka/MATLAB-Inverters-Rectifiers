function x = triangle_generator(N,time,f_sw)
T_sw = 1/f_sw;
x = 0; % initialize wave

% create a 0.5 triangle wave
for n = 1:N*2
    x = (2/(n*pi)^2)*(cos(n*pi)-1)*cos(2*n*pi*f_sw*(time + (T_sw/4))) + x;
end

x = x + 0.5; % add offset of 0.5 so it spans between 0 and 1
