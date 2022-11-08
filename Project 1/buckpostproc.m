function buckpostproc

L = 57.024e-6; %Inductor
C = 473e-6; %Capacitor
R = 1.152; % Rload - LIGHT
D = 1/10; % Duty Cycle

k = 0;

i = 1;
delta_t = (1/10000)/100; %Delta T

initial_I = 3.785; % Initial Inductor Amps
initial_V = 48 + 2; % Initial Capacitor Voltage 
tend = 20*1/10000; %Number of cycles to plot for

[iL, Vload] = buck(L, C, R, delta_t, initial_I, initial_V, tend); %Returns IL and Vload for circuit

while(k(i) < tend)
    sw(i+1) = switching(D, k(i));
    vl(i+1) = (480 - Vload(i+1)).*(sw(i+1)) + (1 - sw(i+1)).*(-Vload(i+1));
    vsw(i+1) = (480)*(1 - sw(i+1));
    id(i+1) = (1 - sw(i+1)).*iL(i+1);
    vd(i+1) = (-480.*sw(i+1));
    isw(i+1) = sw(i+1).*iL(i+1);
    k(i+1) = k(i) + (delta_t);
    i = i + 1;
end

iL_av = avrg(iL, delta_t*100, delta_t);
ic = iL - iL_av;
Pout = (iL - ic) .* Vload;
Pout_av = avrg(Pout, delta_t*100, delta_t);

Pin = 480*isw;
Pin_av = avrg(Pin, delta_t*100, delta_t);
efficency = Pout_av / Pin_av %%Efficiency = 0.9160


%%%% Question 1 Part C, Part IV Code for LIGHT Load

figure(1)
title("Question C, Part IV, Currents (A)")

subplot(5,1,1);
plot(k, (iL - ic), 'r');
ylabel("Output i_R (A)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph
hold on

subplot(5,1,2);
plot(k, ic, 'r');
ylabel("Capacitor i_C (A)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

subplot(5,1,3);
plot(k, isw, 'b');
ylabel("Transistor i_SW (A)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph


subplot(5,1,4);
plot(k, iL, 'c');
ylabel("Inductor i_L (A)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

subplot(5,1,5);
plot(k, id, 'c');
ylabel("Diode i_D (A)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

grid on
hold off

figure(2)
title("Question C, Part IV, Voltages (V)")

subplot(5,1,1);
plot(k, Vload, 'r');
ylabel("Output V_R (V)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph
hold on

subplot(5,1,2);
plot(k, vl, 'r');
ylabel("Inductor V_L (V)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

subplot(5,1,3);
plot(k, vsw, 'b');
ylabel("Transistor V_SW (V)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

subplot(5,1,4);
plot(k, vd, 'c');
ylabel("Diode V_D (V)"); % Y-Label for the graph
xlabel("Time (s)"); % X-Label for the graph

end
