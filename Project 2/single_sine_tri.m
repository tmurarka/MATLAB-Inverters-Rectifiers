% Single Phase Inverter - Sine Triangle

function fdmtl = single_sine_tri(m)

V_dc = 200; % Input DC Voltage
L = 1e-3; % Load Inductance Value
r = 0.1; % Load Resistance Value
f_ac = 200; % Output AC Voltage Frequency
T_ac = 1/f_ac; % Output AC Voltage Period
f_sw = 4000; % Switching Frequency
T_sw = 1/f_sw; % Switching Period 
del_t = T_sw/100; % Time Step
t_end = 10*T_ac; % Simulation End Time

% Initializations
i_ac(1) = -248.7060;
t(1) = 0;
d(1) = 0.5 + 0.5*m*cos(2*pi*200*t(1) - pi/2);
c(1) = d(1) > tri_gen(100,t(1),f_sw);
pos(1) = tri_gen(100,t(1),f_ac) > 0.5;
neg(1) = tri_gen(100,t(1),f_ac) <= 0.5;
V_ag(1) = c(1)*V_dc;
V_bg(1) = (1-c(1))*V_dc;
V_ac(1) = V_ag(1) - V_bg(1);
k = 1;

% Backward Euler Integration Routine
while t(k) < t_end
    d(k+1) = 0.5 + 0.5*m*cos(2*pi*f_ac*(t(k)+del_t) - pi/2);
    c(k+1) = d(k+1) > tri_gen(100,t(k)+del_t,f_sw);
    pos(k+1) = tri_gen(100,t(k)+del_t,f_ac) > 0.5;
    neg(k+1) = tri_gen(100,t(k)+del_t,f_ac) <= 0.5;
    V_ag(k+1) = c(k+1)*V_dc;
    V_bg(k+1) = (1-c(k+1))*V_dc;
    V_ac(k+1) = V_ag(k+1) - V_bg(k+1);
    i_ac(k+1) = (1/(1+(r*del_t/L))) * (i_ac(k) + del_t*V_ac(k+1)/L);
    t(k+1) = t(k) + del_t;
    k = k+1;
end

N = 500;
[avg,ak,bk,rw,err] = fourier(t,V_ac,T_ac,N);

fdmtl = bk(1);

plot(t, V_ac)