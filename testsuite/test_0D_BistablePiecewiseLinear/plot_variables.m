close all
clear all


data = dlmread('results.txt', '', 1, 0);

t = data(:,1);
v = data(:,2);

figure(1)
plot(t,v)
title('Action Potential')
xlabel('Time [ms]')
ylabel('Potential [mV]')


