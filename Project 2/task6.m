% Single Phase Inverter - Sine Triangle
clear all
clc
close all

Pout = [1000];

V_dc = (1/0.95) * sqrt(8000) * sqrt(2); %  input DC
m = (sqrt(8*Pout)*sqrt(2))/V_dc;
L = 1e-3; % inductor
R = 8; % resistor
C = 2e-6; %capacitor
f_ac = 440; % Output AC Voltage Frequency
T_ac = 1/f_ac; % Output AC Voltage Period
w_ac = 2*pi/T_ac;
f_sw = 30000; % Switching Frequency
T_sw = 1/f_sw; % Switching Period 
delta_t = T_sw/100; % Time Step
t_end = 25*T_ac; % Simulation End Time
N = 100;
f = 1/T_ac:1/T_ac:N/T_ac;


% Initializations
k = 1;
i_ac(k) = 0;
Vout(k) = 0; 
t(k) = 0;
d(k) = 0.5 + 0.5*m*sin(w_ac*t(k));
c(k) = d(k) > triangle_generator(N,t(k),f_sw);
V_ag(k) = c(k)*V_dc;
V_bg(k) = (1-c(k))*V_dc;
V_ac(k) = V_ag(k) - V_bg(k);
positive(k) = triangle_generator(N,t(k),f_ac) > 0.5;
negative(k) = triangle_generator(N,t(k),f_ac) <= 0.5;
V_T12(1) = -negative(k)*V_ac(k);
V_T34(1) = positive(k)*V_ac(k);
V_d12(1) = -V_T12(k);
V_d34(1) = -V_T34(k);
i_S12(1) = positive(k)*i_ac(k);
i_S34(1) = -negative(k)*i_ac(k);
i_T12(1) = i_S12(k)*(i_S12(k) > 0);
i_T34(1) = i_S34(k)*(i_S34(k) > 0);
i_d12(1) = i_S12(k)*-(i_S12(k) <= 0);
i_d34(1) = i_S34(k)*-(i_S34(k) <= 0);

% Backward Euler Integration Routine
while t(k) < t_end
    d(k+1) = 0.5 + 0.5*m*sin(w_ac*(t(k)+delta_t));
    c(k+1) = d(k+1) > triangle_generator(N,t(k)+delta_t,f_sw);
    V_ag(k+1) = c(k+1)*V_dc;
    V_bg(k+1) = (1-c(k+1))*V_dc;
    V_ac(k+1) = V_ag(k+1) - V_bg(k+1);
    A_inverse = [1 (delta_t/(L)); (-delta_t/(C)) (1+(delta_t/(R*C)))]^-1;
    D = [delta_t/L; 0]*V_ac(k+1);
    p = A_inverse * ([i_ac(k); Vout(k)] + D);
    i_ac(k+1) = p(1);
    Vout(k+1) = p(2);
    positive(k+1) = triangle_generator(N,t(k)+delta_t,f_ac) > 0.5;
    negative(k+1) = triangle_generator(N,t(k)+delta_t,f_ac) <= 0.5;
    V_T12(k+1) = -negative(k+1)*V_ac(k+1);
    V_T34(k+1) = positive(k+1)*V_ac(k+1);
    V_d12(k+1) = -V_T12(k+1);
    V_d34(k+1) = -V_T34(k+1);
    i_S12(k+1) = positive(k+1)*i_ac(k+1);
    i_S34(k+1) = -negative(k+1)*i_ac(k+1);
    i_T12(k+1) = i_S12(k+1)*(i_S12(k+1) > 0);
    i_T34(k+1) = i_S34(k+1)*(i_S34(k+1) > 0);
    i_d12(k+1) = i_S12(k+1)*-(i_S12(k+1) <= 0);
    i_d34(k+1) = i_S34(k+1)*-(i_S34(k+1) <= 0);
    t(k+1) = t(k) + delta_t;
    k = k+1;
end

[avg,ak,bk,rw,err] = fourser(t,Vout,T_ac,N);

thd = THD(ak, bk, N)

for k = 1:N
    ck(k) = ((ak(k))^2 + (bk(k))^2)^(1/2);
end

fundamental = ck(1)

figure;
subplot(2,1,1)
plot(t,V_ac)
xlim([23*T_ac 24*T_ac])
ylim([-300 300])
title("AC Voltage")
xlabel("t (s)")
ylabel("V_a_c (V)")

subplot(2,1,2)
plot(t,Vout)
xlim([23*T_ac 24*T_ac])
ylim([-300 300])
title("Output Voltage")
xlabel("t (s)")
ylabel("V_o_u_t (V)")

figure;
plot(t,i_ac)
xlim([23*T_ac 24*T_ac])
ylim([-25 25])
title("Output AC Current")
xlabel("t (s)")
ylabel("i_a_c (A)")

figure;
stem(f,ck)
title("ck Vs frequency")
xlabel("Frequency (Hz)")
ylabel("ck")

figure;
subplot(4,2,1);
plot(t,V_T12)
xlim([23*T_ac 24*T_ac])
title("Transistor 1 = 2 Voltage")
xlabel("t (s)")
ylabel("V_T_1_=_T_2 (V)")

subplot(4,2,2);
plot(t,V_T34)
xlim([23*T_ac 24*T_ac])
title("Transistor 3 = 4 Voltage")
xlabel("t (s)")
ylabel("V_T_3_=_T_4 (V)")

subplot(4,2,3);
plot(t,V_d12)
xlim([23*T_ac 24*T_ac])
title("Diode 1 = 2 Voltage")
xlabel("t (s)")
ylabel("V_d_1_=_d_2 (V)")

subplot(4,2,4);
plot(t,V_d34)
xlim([23*T_ac 24*T_ac])
title("Diode 3 = 4 Voltage")
xlabel("t (s)")
ylabel("V_d_3_=_d_4 (V)")

subplot(4,2,5);
plot(t,i_T12)
xlim([23*T_ac 24*T_ac])
title("Transistor 1 = 2 Current")
xlabel("t (s)")
ylabel("i_T_1_=_T_2 (A)")

subplot(4,2,6);
plot(t,i_T34)
xlim([23*T_ac 24*T_ac])
title("Transistor 3 = 4 Current")
xlabel("t (s)")
ylabel("i_T_3_=_T_4 (A)")

subplot(4,2,7);
plot(t,i_d12)
xlim([23*T_ac 24*T_ac])
title("Diode 1 = 2 Current")
xlabel("t (s)")
ylabel("i_d_1_=_d_2 (A)")

subplot(4,2,8);
plot(t,i_d34)
xlim([23*T_ac 24*T_ac])
title("Diode 3 = 4 Current")
xlabel("t (s)")
ylabel("i_d_3_=_d_4 (A)")