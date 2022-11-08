clear all
clc
close all

% given quantities
L_ac = 2e-3; % ac inductance source
L_dc = 6e-3; % dc inductance
R_load = 1.2471; % load resistance
f_ac = 60; % frequency of fundamental
T_ac = 1/f_ac; % period of fundamental
w_ac = 2*pi/T_ac;
E_given = sqrt(2)*120; % input voltage
delta_t = 10^-7; % Time Step
t_end = 9*T_ac; % Simulation End Time
k = 1;
t(k) = 0;
v_dc(k) = 0;
i_dc(k) = 0;
tau = 10^-5;
epsilon = .1;
i_a(k) = 0;
i_b(k) = 0;
i_c(k) = 0;

while t(k) < t_end
    e_as(k+1) = E_given*cos(w_ac*t(k));
	e_bs(k+1) = E_given*cos((w_ac*t(k)) - ((2*pi)/3));
	e_cs(k+1) = E_given*cos((w_ac*t(k)) + ((2*pi)/3));

	if i_a(k) > 0 
		i_d1(k+1) = i_a(k);
	else 
		i_d1(k+1) = 0;
	end

	if i_b(k) > 0
		i_d3(k+1) = i_b(k);
	else
		i_d3(k+1) = 0;
	end

	if i_c(k) > 0
		i_d5(k+1) = i_c(k);
	else
		i_d5(k+1) = 0;
	end

	i_dc(k+1) = i_d1(k+1) + i_d3(k+1)+ i_d5(k+1);

	v_dc(k+1) = (tau/(tau+delta_t)) * (v_dc(k) + ((L_dc/tau) * (i_dc(k+1) - i_dc(k))) + ((R_load*i_dc(k+1)*delta_t)/tau));

	if i_a(k) >= epsilon
		v_aL(k+1) = v_dc(k+1);
    elseif i_a(k) <= (-1*epsilon)
		v_aL (k+1) = 0;
	else
		v_aL(k+1) = ((v_dc(k+1)*i_a(k))/(2*epsilon)) + (v_dc(k+1)/2);
	end

	if i_b(k) >= epsilon
		v_bL(k+1) = v_dc(k+1);
    elseif i_b(k) <= (-1*epsilon)
		v_bL(k+1) = 0;
	else
		v_bL(k+1) = ((v_dc(k+1)*i_b(k))/(2*epsilon)) + (v_dc(k+1)/2);
	end

	if i_c(k) >= epsilon
		v_cL(k+1) = v_dc(k+1);
    elseif i_c(k) <= (-1*epsilon)
		v_cL(k+1) = 0;
	else
		v_cL(k+1) = ((v_dc(k+1)*i_c(k))/(2*epsilon)) + (v_dc(k+1)/2);
	end 

	va(k+1) = ((2/3)*v_aL(k+1)) - ((1/3)*v_bL(k+1)) - ((1/3)*v_cL(k+1));
	vb(k+1) = ((2/3)*v_bL(k+1)) - ((1/3)*v_cL(k+1)) - ((1/3)*v_aL(k+1));
	vc(k+1) = ((2/3)*v_cL(k+1)) - ((1/3)*v_aL(k+1)) - ((1/3)*v_bL(k+1));

	i_a(k+1) = i_a(k) + ((delta_t*(e_as(k+1)-va(k+1)))/L_ac);
	i_b(k+1) = i_b(k) + ((delta_t*(e_bs(k+1)-vb(k+1)))/L_ac);
	i_c(k+1) = i_c(k) + ((delta_t*(e_cs(k+1)-vc(k+1)))/L_ac);
		
	t(k+1) = t(k) + delta_t;
	
    k = k+1;
end

for N = 1:t_end/T_ac
    v_dcavg(N) = avrg(v_dc,N*T_ac,delta_t);
    i_dcavg(N) = avrg(i_dc,N*T_ac,delta_t);
end

theta = t * w_ac;

figure;
subplot(2,1,1)
plot(theta*180/pi,v_dc)
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
title("DC Voltage (V)")
xlabel("{\theta_a_c} (deg)")
ylabel("V_d_c (V)")

subplot(2,1,2)
plot(theta*180/pi,i_dc)
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
title("DC Current (A)")
xlabel("{\theta_a_c} (deg)")
ylabel("i_d_c (A)")

figure;
subplot(4,1,1)
plot(theta*180/pi, i_a, 'b')
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
xlabel("{\theta_a_c} (deg)")
ylabel("I_a_c (A)")
title("Phase Currents for Q2 - Mode 2")

subplot(4,1,2)
plot(theta*180/pi, i_b, 'r')
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
xlabel("{\theta_a_c} (deg)")
ylabel("I_b_c (A)")

subplot(4,1,3)
plot(theta*180/pi, i_c, 'c')
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
xlabel("{\theta_a_c} (deg)")
ylabel("I_c_c (A)")

subplot(4,1,4)
plot(theta*180/pi, i_a+i_b+i_c, 'g')
xlim([(t_end-T_ac)*w_ac*180/pi t_end*w_ac*180/pi])
xlabel("{\theta_a_c} (deg)")
ylabel("I_t_o_t_a_l (A)")

figure;
plot(i_dcavg,v_dcavg)
xlabel("i_d_c_ _a_v_g (A)")
ylabel("V_d_c_ _a_v_g (V)")
title("V_d_c vs i_d_c for Q2 - Mode 2")
