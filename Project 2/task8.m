% 3 Phase Inverter
clear all
clc
close all

r = 0.5;
V_dc = 80;
L = 5e-3;
f_ac = 60;
T_ac = 1/f_ac;
del_t = T_ac/100; % Time Step
t_end = 25*T_ac; % Simulation End Time
N = 50;
f = 1/T_ac:1/T_ac:N/T_ac;

k = 1;
i_as(k) = 20;
i_bs(k) = 20*cos(-2*pi/3);
i_cs(k) = 20*cos(2*pi/3);

t(k) = 0;
d_a(k) = 0.5 + 0.5*cos(2*pi*f_ac*t(k));
d_b(k) = 0.5 + 0.5*cos(2*pi*f_ac*t(k) - 2*pi/3);
d_c(k) = 0.5 + 0.5*cos(2*pi*f_ac*t(k) + 2*pi/3);

c_a(k) = d_a(k) > tri_gen(N,t(k),f_ac);
c_b(k) = d_b(k) > tri_gen(N,t(k),f_ac);
c_c(k) = d_c(k) > tri_gen(N,t(k),f_ac);

V_ag(k) = c_a(k)*V_dc;
V_bg(k) = c_b(k)*V_dc;
V_cg(k) = c_c(k)*V_dc;

V_as(k) = (2/3)*V_ag(k) - (1/3)*V_bg(k) - (1/3)*V_cg(k);
V_bs(k) = (2/3)*V_bg(k) - (1/3)*V_ag(k) - (1/3)*V_cg(k);
V_cs(k) = (2/3)*V_cg(k) - (1/3)*V_ag(k) - (1/3)*V_bg(k);

while t(k) < t_end 
   
    d_a(k+1) = 0.5 + 0.5*cos(2*pi*f_ac*t(k));
    d_b(k+1) = 0.5 + 0.5*cos(2*pi*f_ac*t(k) - 2*pi/3);
    d_c(k+1) = 0.5 + 0.5*cos(2*pi*f_ac*t(k) + 2*pi/3);

    c_a(k+1) = d_a(k+1) > tri_gen(N,t(k)+del_t,f_ac);
    c_b(k+1) = d_b(k+1) > tri_gen(N,t(k)+del_t,f_ac);
    c_c(k+1) = d_c(k+1) > tri_gen(N,t(k)+del_t,f_ac);

    V_ag(k+1) = c_a(k+1)*V_dc;
    V_bg(k+1) = c_b(k+1)*V_dc;
    V_cg(k+1) = c_c(k+1)*V_dc;

    V_as(k+1) = (2/3)*V_ag(k+1) - (1/3)*V_bg(k+1) - (1/3)*V_cg(k+1);
    V_bs(k+1) = (2/3)*V_bg(k+1) - (1/3)*V_ag(k+1) - (1/3)*V_cg(k+1);
    V_cs(k+1) = (2/3)*V_cg(k+1) - (1/3)*V_ag(k+1) - (1/3)*V_bg(k+1);
    
    i_as(k+1) = i_as(k) + del_t*(V_as(k+1) - r*i_as(k))/L;
    i_bs(k+1) = i_bs(k) + del_t*(V_bs(k+1) - r*i_bs(k))/L;
    i_cs(k+1) = i_cs(k) + del_t*(V_cs(k+1) - r*i_cs(k))/L;
    
    t(k+1) = t(k) + del_t;
    k = k + 1;
end

figure;
subplot(3,1,1)
plot(t,i_as)
hold on
xlim([4*T_ac 5*T_ac])
title("Phase a current")
xlabel("t (s)")
ylabel("i_a_s (A)")

subplot(3,1,2)
plot(t,i_bs)
hold on
xlim([4*T_ac 5*T_ac])
title("Phase b current")
xlabel("t (s)")
ylabel("i_b_s (A)")

subplot(3,1,3)
plot(t,i_cs)
hold on
xlim([4*T_ac 5*T_ac])
title("Phase c current")
xlabel("t (s)")
ylabel("i_c_s (A)")

figure;
subplot(3,1,1)
plot(t,V_ag)
xlim([4*T_ac 5*T_ac])
title("V_a_g Voltage")
xlabel("t (s)")
ylabel("V_a_s (V)")

subplot(3,1,2)
plot(t,V_bg)
xlim([4*T_ac 5*T_ac])
title("V_b_g Voltage")
xlabel("t (s)")
ylabel("V_b_s (V)")

subplot(3,1,3)
plot(t,V_cg)
xlim([4*T_ac 5*T_ac])
title("V_c_g Voltage")
xlabel("t (s)")
ylabel("V_c_s (V)")

figure;
subplot(3,1,1)
plot(t,V_as)
xlim([4*T_ac 5*T_ac])
title("V_a_s Voltage")
xlabel("t (s)")
ylabel("V_a_s (V)")

subplot(3,1,2)
plot(t,V_bs)
xlim([4*T_ac 5*T_ac])
title("V_b_s Voltage")
xlabel("t (s)")
ylabel("V_b_s (V)")

subplot(3,1,3)
plot(t,V_cs)
xlim([4*T_ac 5*T_ac])
title("V_c_s Voltage")
xlabel("t (s)")
ylabel("V_c_s (V)")
