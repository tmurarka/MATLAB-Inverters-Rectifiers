function [av] = avrg(x,T,dt)
    n_period = T/dt;
    av = 0;

for m = 1:n_period
    av = av + x(m)*dt/T;
end

