clear all
clc
close all

T = 2*pi/(2*pi*100); % fundamental frequency
t = 0:T/1000:1; % time
x = 10*cos(2*pi*100*t) + 30*sin(2*pi*300*t)-5*cos(2*pi*1000*t)+17*sin(2*pi*1500*t)-5*cos(2*pi*5000*t); % input
N = 50; % terms

[avg,ak,bk,rw,err] = fourser(t,x,T,N);

harmonic = THD(ak, bk, N)

figure;
subplot(2,1,1);
plot(t,rw)
xlim([3*T 4*T])
title("Fourier Output")
xlabel("t(s)")
ylabel("Fourier")

subplot(2,1,2);
plot(t,x)
xlim([3*T 4*T])
title("Original")
xlabel("t(s)")
ylabel("Original")

err;
