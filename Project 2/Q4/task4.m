clear all
clc
close all

Vac_given = 400/pi; % V_ac_given
Vdc = Vac_given*pi/4; % Vdc

% given quantities
R = 0.5; % resistor
L = 1e-3; % inductor
f_ac = 60; % frequency of fundamental
T_ac = 1/f_ac; % period of fundamental
N = 30; %number of terms for fourier
frequency = 1/T_ac:1/T_ac:N/T_ac; % creates the frequency space for ck values

tau = L/R; % time constant

I_max = (Vdc/R)*(1-exp(-T_ac/(2*tau)))/(1+exp(-T_ac/(2*tau))); % maximum ac
I_min = -I_max; % minimum ac

t = 0:T_ac/1000:T_ac; % time

for k = 1:length(t) % iteration
    if t(k) <= T_ac/2 % first half of interval
        V_ac(k) = Vdc; 
        V_s12(k) = 0;
        V_s34(k) = V_ac(k);
        i_ac(k) = Vdc/R + (I_min - (Vdc/R))*exp(-t(k)/tau);
        i_s12(k) = i_ac(k);
        i_s34(k) = 0;
        i_dc(k) = i_ac(k);
    else
        V_ac(k) = -Vdc;
        V_s12(k) = -V_ac(k);
        V_s34(k) = 0;
        i_ac(k) = -Vdc/R + (I_max + (Vdc/R))*exp(-(t(k)-(T_ac/2))/tau);
        i_s12(k) = 0;
        i_s34(k) = -i_ac(k);
        i_dc(k) = -i_ac(k);
    end
end

[avg,ak,bk,rw,err] = fourser(t,V_ac,T_ac,N); % fits a fourier function to V_ac

harmonic = THD(ak, bk, N) % finds harmonics of V_ac

% Transistor Currents
i_T12 = i_s12 .* (i_s12 >= 0); % ensures S12 current is greater than 0
i_T34 = i_s34 .* (i_s34 >= 0); % ensures S34 current is greater than 0

% Diode Currents
i_d12 = i_s12 .* -(i_s12 <= 0);
i_d34 = i_s34 .* -(i_s34 <= 0);

figure;
subplot(3,2,1);
plot(t,i_dc)
ylim([-250 250])
title("DC Current")
xlabel("t (s)")
ylabel("i_d_c (A)")

subplot(3,2,2);
plot(t,i_T12)
ylim([-50 250])
title("Transistor 1 = 2 Current")
xlabel("t (s)")
ylabel("i_T_1_=_T_2 (A)")

subplot(3,2,3);
plot(t,i_d12)
ylim([-50 250])
title("Diode 1 = 2 Current")
xlabel("t (s)")
ylabel("i_d_1_=_d_2 (A)")

subplot(3,2,4);
plot(t,i_d34)
ylim([-50 250])
title("Diode 3 = 4 Current")
xlabel("t (s)")
ylabel("i_d_3_=_d_4 (A)")

subplot(3,2,5);
plot(t,i_T34)
ylim([-50 250])
title("Transistor 3 = 4 Current")
ylabel("i_T_3_=_T_4 (A)")
xlabel("t (s)")

subplot(3,2,6);
plot(t,i_ac)
ylim([-250 250])
title("AC Current")
xlabel("t(s)")
ylabel("i_a_c (A)")
