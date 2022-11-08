clear all
clc

% given quantities
L_ac = 2e-3; % ac inductor
max_gamma = pi/3;
max_alpha = pi/6;
max_delta = pi/3;
f_ac = 60; % frequency of fundamental
T_ac = 1/f_ac; % period of fundamental
w_ac = 2*pi/T_ac;
gamma = 0:max_gamma/1000:max_gamma;
alpha = 0:max_alpha/1000:max_alpha;
delta = 0:max_delta/1000:max_delta;
E_given = sqrt(2)*120;

%% Part A
k = 0;
for N = 1:length(gamma)
    if (gamma(N) >= 0 && gamma(N) <= pi/3)
        V_dc_1(k+1) = (3*sqrt(3)*E_given/(2*pi))*(1+cos(gamma(N)));
        i_dc_1(k+1) = -(30000*cos(gamma(N))*sqrt(6)/377)+(30000*sqrt(6)/377);
    k = k+1;
    end
end
k = 0;
figure;
plot(i_dc_1,V_dc_1, 'b')
hold on

for N = 1:length(alpha)
    if (alpha(N) >= 0 && alpha(N) <= pi/6)
        V_dc_2(k+1) = 9*E_given*cos(alpha(N) + (pi/6))/(2*pi);
        i_dc_2(k+1) = sqrt(3)*E_given*sin(alpha(N) + (pi/6))/(2*w_ac*L_ac);
    k = k+1;
    end
end
k = 0;

plot(i_dc_2,V_dc_2, 'r')

for N = 1:length(delta)
    if (delta(N) >= 0 && delta(N) <= pi/3)
        V_dc_3(k+1) = 9*E_given*(1-sin(delta(N) + pi/6))/(2*pi);
        i_dc_3(k+1) = (1+sin(delta(N) + pi/6))*E_given/(2*w_ac*L_ac);
    k = k+1;
    end
end

plot(i_dc_3,V_dc_3, 'g')
ylim([-25 300])
xlim([-10 250])
xlabel("i_d_c (A)")
ylabel("V_d_c (V)")
title("V_d_c vs i_d_c for Q1-A")

%% Part B

k = 0;
max_gamma = pi/12;
max_theta = 2*pi;
theta = 0:max_theta/1000:max_theta;

V_dc = (3*sqrt(3)*E_given/(2*pi))*(1+cos(pi/12)); % 275.9087 V
i_dc = -(30000*cos(pi/12)*sqrt(6)/377)+(30000*sqrt(6)/377); % 6.6417 A
R_load_PartB = V_dc/i_dc; % 41.5417 Ohms

for N = 1:length(theta)
    if (theta(N) >= 0 && theta(N) <= max_gamma) % Commutation
        i_ac(k+1) = i_dc;
        i_bc(k+1) = ((sqrt(3)*(1-cos(theta(N)))*E_given)/(2*w_ac*L_ac))-i_dc;
        i_cc(k+1) = -i_bc(k+1) - i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_gamma && theta(N) <= pi/3)
        i_ac(k+1) = i_dc;
        i_bc(k+1) = 0;
        i_cc(k+1) = -i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > pi/3 && theta(N) <= max_gamma + pi/3) % Commutation
        i_ac(k+1) = (E_given*(-sin((3*theta(N)-2*pi)/3)+sin(theta(N))-sqrt(3)))/(2*L_ac*w_ac) + i_dc;
        i_bc(k+1) = -i_ac(k+1) + i_dc;
        i_cc(k+1) = -i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_gamma + pi/3 && theta(N) <= 2*pi/3)
        i_ac(k+1) = 0;
        i_bc(k+1) = i_dc;
        i_cc(k+1) = -i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > 2*pi/3 && theta(N) <= max_gamma + 2*pi/3) % Commutation
        i_ac(k+1) = (E_given*(-sin((3*theta(N)+2*pi)/3)+sin(theta(N))-sqrt(3)))/(2*L_ac*w_ac);
        i_bc(k+1) = i_dc;
        i_cc(k+1) = -i_ac(k+1) - i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_gamma + 2*pi/3 && theta(N) <= 3*pi/3) 
        i_ac(k+1) = -i_dc;
        i_bc(k+1) = i_dc;
        i_cc(k+1) = 0;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > 3*pi/3 && theta(N) <= max_gamma + 3*pi/3) % Commutation
        i_ac(k+1) = -i_dc;
        i_bc(k+1) = (sqrt(3)*E_given*(-cos(theta(N))-1))/(2*L_ac*w_ac) + i_dc;
        i_cc(k+1) = -i_bc(k+1) + i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_gamma + 3*pi/3 && theta(N) <= 4*pi/3) 
        i_ac(k+1) = -i_dc;
        i_bc(k+1) = 0;
        i_cc(k+1) = i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > 4*pi/3 && theta(N) <= max_gamma + 4*pi/3) % Commutation
        i_bc(k+1) = (E_given*(sin((3*theta(N)-2*pi)/3)-sin(theta(N))-sqrt(3)))/(2*L_ac*w_ac);
        i_ac(k+1) = -i_bc(k+1) - i_dc;
        i_cc(k+1) = i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_gamma + 4*pi/3 && theta(N) <= 5*pi/3) 
        i_ac(k+1) = 0;
        i_bc(k+1) = -i_dc;
        i_cc(k+1) = i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > 5*pi/3 && theta(N) <= max_gamma + 5*pi/3) % Commutation
        i_bc(k+1) = -i_dc;
        i_cc(k+1) = (E_given*(sin((3*theta(N)+2*pi)/3)-sin(theta(N))-sqrt(3)))/(2*L_ac*w_ac) + i_dc;
        i_ac(k+1) = -i_cc(k+1) + i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
 
    if(theta(N) > max_gamma + 5*pi/3 && theta(N) <= 2*pi)
        i_ac(k+1) = i_dc;
        i_bc(k+1) = -i_dc;
        i_cc(k+1) = 0;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    k = k+1;
end

figure;
subplot(4,1,1)
plot(theta*180/pi, i_ac, 'b')
xlabel("{\theta_a_c}")
ylabel("I_a_c (A)")
title("Phase Currents for Q1-B")

subplot(4,1,2)
plot(theta*180/pi, i_bc, 'r')
xlabel("{\theta_a_c}")
ylabel("I_b_c (A)")

subplot(4,1,3)
plot(theta*180/pi, i_cc, 'c')
xlabel("{\theta_a_c}")
ylabel("I_c_c (A)")

subplot(4,1,4)
plot(theta*180/pi, i_total, 'g')
xlabel("{\theta_a_c}")
ylabel("I_t_o_t_a_l (A)")



%% Part C

k = 0;
max_alpha = pi/12;
max_theta = max_alpha + 6*pi/3;
theta = max_alpha:max_alpha/1000:max_theta;

V_dc = (9*E_given/(2*pi))*cos(max_alpha+(pi/6)); % 171.8873 V
i_dc = (sqrt(3)*E_given/(2*w_ac*L_ac))*sin(max_alpha+(pi/6)); % 137.8322 A
R_load_PartC = V_dc/i_dc; % 1.2471 Ohms

for N = 1:length(theta)
    if(theta(N) >= max_alpha && theta(N) <= max_alpha + pi/3) % Commutation
        i_ac(k+1) = i_dc;
        i_bc(k+1) = (sqrt(3)*E_given*(cos(pi/12)-cos(theta(N))))/(2*L_ac*w_ac)-i_dc;
        i_cc(k+1) = -i_bc(k+1) - i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_alpha + pi/3 && theta(N) <= max_alpha + 2*pi/3) % Commutation
        i_ac(k+1) = -(E_given*(sqrt(2)*(sin((3*theta(N)-2*pi)/3)-sin(theta(N))+sin((5*pi)/12))+1))/(2^(3/2)*L_ac*w_ac) + i_dc;
        i_bc(k+1) = -i_ac(k+1) + i_dc;
        i_cc(k+1) = -i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_alpha + 2*pi/3 && theta(N) <= max_alpha + 3*pi/3) % Commutation
        i_ac(k+1) = -(E_given*(sqrt(2)*(sin((3*theta(N)+2*pi)/3)-sin(theta(N))-sin((17*pi)/12))+1))/(2^(3/2)*L_ac*w_ac);
        i_bc(k+1) = i_dc;
        i_cc(k+1) = -i_ac(k+1) - i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_alpha + 3*pi/3 && theta(N) <= max_alpha + 4*pi/3) % Commutation
        
        i_ac(k+1) = -i_dc;
        i_bc(k+1) = (sqrt(3)*E_given*(cos((13*pi)/12)-cos(theta(N))))/(2*L_ac*w_ac) + i_dc;
        i_cc(k+1) = -i_bc(k+1) + i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_alpha + 4*pi/3 && theta(N) <= max_alpha + 5*pi/3) % Commutation
        i_bc(k+1) = (E_given*(sqrt(2)*(sin((3*theta(N)-2*pi)/3)-sin(theta(N))+sin((17*pi)/12))-1))/(2^(3/2)*L_ac*w_ac);
        i_ac(k+1) = -i_bc(k+1) - i_dc;
        i_cc(k+1) = i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    if(theta(N) > max_alpha + 5*pi/3 && theta(N) <= max_alpha + 6*pi/3) % Commutation
        i_bc(k+1) = -i_dc;
        i_cc(k+1) = (E_given*(sqrt(2)*(sin((3*theta(N)+2*pi)/3)-sin(theta(N))-sin((29*pi)/12))-1))/(2^(3/2)*L_ac*w_ac) + i_dc;
        i_ac(k+1) = -i_cc(k+1) + i_dc;
        i_total(k+1) = i_ac(k+1) + i_bc(k+1) + i_cc(k+1);
    end
    
    k = k+1;
end


figure;
subplot(4,1,1)
plot(theta*180/pi, i_ac, 'b')
xlabel("{\theta_a_c} (deg)")
ylabel("I_a_c (A)")
title("Phase Currents for Q1-C")

subplot(4,1,2)
plot(theta*180/pi, i_bc, 'r')
xlabel("{\theta_a_c} (deg)")
ylabel("I_b_c (A)")

subplot(4,1,3)
plot(theta*180/pi, i_cc, 'c')
xlabel("{\theta_a_c} (deg)")
ylabel("I_c_c (A)")

subplot(4,1,4)
plot(theta*180/pi, i_total, 'g')
xlabel("{\theta_a_c} (deg)")
ylabel("I_t_o_t_a_l (A)")
