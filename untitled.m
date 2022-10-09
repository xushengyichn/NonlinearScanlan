clc
clear
close all
t=0:0.01:10;
f=1;
omega=2*pi*f;

y=sin(omega*t);
figure
plot(t,y);

L=0:0.01:1;

phi=sin(pi*L/max(L));
figure
plot(L,phi)

phi'.*y