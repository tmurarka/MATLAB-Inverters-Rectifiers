clear all
clc
close all

fundamental = [38.4228 54.2272 65.9464 76.3693 85.3511 120.8963]
m = [0.3004    0.4249    0.5203    0.6008    0.6718    0.9500]

figure;
plot(m, fundamental)
xlim([0 1])
ylim([0 150])
title("fundamental vs varying duty cycle")
xlabel("m (duty cycle)")
ylabel("fundamental (ck)")
