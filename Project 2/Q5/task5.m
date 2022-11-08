% Single Phase Inverter - 180 degree Switching

clear
clc
close all

V_dc = 100; %DC Voltage
L = 1e-3; % inductor
r = 0.5; % resistor
f_ac = 60; % switching frequency
T_ac = 1/f_ac; % switching period
delta_t = T_ac/1000; % time stepping
t_end = 25*T_ac; % end time
N = 35; % number of terms
frequency = 1/T_ac:1/T_ac:N/T_ac; %frequency space
k = 1;


i_ac(k) = -193.8931384478502; % initial I_ac from part 4
t(k) = 0;
T12(k) = 1; % assume T12 are on
T34(k) = 0; % assume T23 are off
V_ag(k) = T12(k)*V_dc; %calculates Va phase leg
V_bg(k) = T34(k)*V_dc; % calcualtes Vb phase leg
V_ac(k) = V_ag(k) - V_bg(k); % Vac
V_T12(k) = -T34(k)*V_ac(k); % Transistor voltage
V_T34(k) = T12(k)*V_ac(k);
V_d12(k) = -V_T12(k);
V_d34(k) = -V_T34(k);
i_S12(k) = T12(k)*i_ac(k); % switch current
i_S34(k) = -T34(k)*i_ac(k);
i_T12(k) = i_S12(k)*(i_S12(k) > 0); % switch current is transistor current
i_T34(k) = i_S34(k)*(i_S34(k) > 0);
i_d12(k) = i_S12(k)*-(i_S12(k) <= 0); % switch current is diode current
i_d34(k) = i_S34(k)*-(i_S34(k) <= 0);
i_dc(k) = i_S12(k) + i_S34(k);

% Backward Euler
while t(k) < t_end
    T12(k+1) = triangle_generator(N,t(k)+delta_t,f_ac) > 0.5; % generates triangle wave to be compared for half interval
    T34(k+1) = triangle_generator(N,t(k)+delta_t,f_ac) <= 0.5;
    V_ag(k+1) = T12(k+1)*V_dc;
    V_bg(k+1) = T34(k+1)*V_dc;
    V_ac(k+1) = V_ag(k+1) - V_bg(k+1);
    i_ac(k+1) = (1/(1+(r*delta_t/L))) * (i_ac(k) + delta_t*V_ac(k+1)/L);
    V_T12(k+1) = -T34(k+1)*V_ac(k+1);
    V_T34(k+1) = T12(k+1)*V_ac(k+1);
    V_d12(k+1) = -V_T12(k+1);
    V_d34(k+1) = -V_T34(k+1);
    i_S12(k+1) = T12(k+1)*i_ac(k+1);
    i_S34(k+1) = -T34(k+1)*i_ac(k+1);
    i_T12(k+1) = i_S12(k+1)*(i_S12(k+1) > 0);
    i_T34(k+1) = i_S34(k+1)*(i_S34(k+1) > 0);
    i_d12(k+1) = i_S12(k+1)*-(i_S12(k+1) <= 0);
    i_d34(k+1) = i_S34(k+1)*-(i_S34(k+1) <= 0);
    i_dc(k+1) = i_S12(k+1) + i_S34(k+1);
    t(k+1) = t(k) + delta_t;
    k = k+1;
end

[avg,ak,bk,rw,err] = fourser(t,V_ac,T_ac,N);

thd = THD(ak, bk, N)

for k = 1:N
    ck(k) = ((ak(k))^2 + (bk(k))^2)^(1/2);
end


figure;
subplot(4,2,1);
plot(t,V_T12)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Transistor 1 and 2 Voltage")
xlabel("t (s)")
ylabel("V_T_1_=_T_2 (V)")

subplot(4,2,2);
plot(t,i_T12)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Transistor 1 and 2 Current")
xlabel("t (s)")
ylabel("i_T_1_=_T_2 (A)")

subplot(4,2,3);
plot(t,V_T34)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Transistor 3 and 4 Voltage")
xlabel("t (s)")
ylabel("V_T_3_=_T_4 (V)")

subplot(4,2,4);
plot(t,i_T34)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Transistor 3 and 4 Current")
xlabel("t (s)")
ylabel("i_T_3_=_T_4 (A)")


subplot(4,2,5);
plot(t,V_d12)
xlim([20*T_ac 21*T_ac])
ylim([-250 50])
title("Diode 1 and 2 Voltage")
xlabel("t (s)")
ylabel("V_d_1_=_d_2 (V)")

subplot(4,2,6);
plot(t,i_d12)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Diode 1 and 2 Current")
xlabel("t (s)")
ylabel("i_d_1_=_d_2 (A)")


subplot(4,2,7);
plot(t,V_d34)
xlim([20*T_ac 21*T_ac])
ylim([-250 50])
title("Diode 3 and 4 Voltage")
xlabel("t (s)")
ylabel("V_d_3_=_d_4 (V)")

subplot(4,2,8);
plot(t,i_d34)
xlim([20*T_ac 21*T_ac])
ylim([-50 250])
title("Diode 3 and 4 Current")
xlabel("t (s)")
ylabel("i_d_3_=_d_4 (A)")

figure;
stem(frequency,ck)
title("ck Vs frequency")
xlabel("Frequency (Hz)")
ylabel("ck")

figure;
subplot(2,1,1);
plot(t,rw)
xlim([20*T_ac 22*T_ac])
title("Reconstructed Output")
xlabel("t (s)")
ylim([-150 150])
ylabel("V_a_c(Reconstructed) (V)")

subplot(2,1,2);
plot(t,V_ac)
xlim([20*T_ac 22*T_ac])
ylim([-150 150])
title("Original Output Voltage")
xlabel("t (s)")
ylabel("V_a_c(Original) (V)")

figure;
subplot(2,1,1)
plot(t,i_ac)
xlim([20*T_ac 21*T_ac])
ylim([-250 250])
title("Output AC Current")
xlabel("t (s)")
ylabel("i_a_c (A)")

subplot(2,1,2)
plot(t,i_dc)
xlim([20*T_ac 21*T_ac])
ylim([-250 250])
title("DC Current")
xlabel("t (s)")
ylabel("i_d_c (A)")