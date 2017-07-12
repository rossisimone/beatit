close all
clear all


data = dlmread('output.txt', '', 1, 0);

t = data(:,1);
v = data(:,2);
%cai = data(:,7);

figure(1)
plot(t,v,'-r')
title('Action Potential')
xlabel('Time [ms]')
ylabel('Potential [mV]')


data2 = dlmread('results.txt', '', 1, 0);

t2 = data(:,1);
v2 = data(:,2);
%cai = data(:,7);

figure(1)
hold on
plot(t2,v2, '-b')
xlim([0, 600])


max(abs(v-v2))
