close all
clear all


data = dlmread('output.txt', '', 1, 0);

t = data(:,1);
v = data(:,2);
cai = data(:,7);

figure(1)
plot(t,v)
title('Action Potential')
xlabel('Time [ms]')
ylabel('Potential [mV]')


figure(2)
plot(t,cai)
title('Calcium Concentration')
xlabel('Time [ms]')
ylabel('Cai [mM]')



data = dlmread('results.txt', '', 1, 0);

t = data(:,1);
v = data(:,2);
cai = data(:,7);

figure(1)
plot(t,v)
title('Action Potential')
xlabel('Time [ms]')
ylabel('Potential [mV]')


figure(2)
plot(t,cai)
title('Calcium Concentration')
xlabel('Time [ms]')
ylabel('Cai [mM]')
