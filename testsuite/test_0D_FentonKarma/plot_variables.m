close all
clear all


data = dlmread('results.txt', '', 1, 0);

t = data(:,1);
V = data(:,2);
v = data(:,3);
w = data(:,4);

figure(1)
plot(t,V)
title('Action Potential')
xlabel('Time [ms]')
ylabel('Potential [mV]')

figure(2)
plot(t,v)
title('Recovery Variable v')
xlabel('Time [ms]')
ylabel('v')

figure(3)
plot(t,w)
title('Recovery Variable w')
xlabel('Time [ms]')
ylabel('w')

