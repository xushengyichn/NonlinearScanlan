clc
clear
close all

load doubleTMDcontrol_bad.mat
% load singleTMDcontrol.mat
t=0:0.01:t_length;


figure
subplot(3,1,1)
plot(t,dis_beam)
title('Displacement of Beam')
xlabel('Time (s)')
ylabel('Displacement (m)')
subplot(3,1,2)
plot(t,dis_TMD1_all)
title('Displacement of TMD1')
xlabel('Time (s)')
ylabel('Displacement (m)')
subplot(3,1,3)
plot(t,dis_TMD2_all)
title('Displacement of TMD2')
xlabel('Time (s)')
ylabel('Displacement (m)')

figure

fs = 100;
t = 0:1/fs:t_length;
x = dis_beam;

pspectrum(x,fs,'spectrogram', ...
    'TimeResolution',60,'Overlap',86,'Leakage',0.875)

dis_beam=dis_beam(end*3/4:end*4/4)
[psd_avg, f, psd_plot] =fft_transfer(fs,dis_beam')
figure
plot(f,psd_plot)