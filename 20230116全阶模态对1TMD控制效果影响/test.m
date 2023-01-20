
data=u(1,:);
t=0.01;
fs=1/t;

[psd_avg, f, psd_plot] = fft_transfer_modified(fs,data');
figure
plot(f,psd_plot)