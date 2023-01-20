function [psd_avg, f, psd_plot] = fft_transfer_modified(fs,data)
M = 1;
N = fix(size(data,1)/M);
f = fs*(1:fix(N/2))/N;
f = f';
psd_out = cell(N,M);
for k1 = 1:M
    Y = fft(data(((k1-1)*N+1):k1*N,:));
    for k2 = 2:length(f)+1
        psd_temp = sqrt(1/fs/(2*pi*N))*Y(k2,:);
        psd_out{k2-1,k1} = transpose(psd_temp)*conj(psd_temp);
    end
    %     disp((k1)/M);
end

psd_avg = cell(length(f),1);
for k2 = 1:length(f)
    psd_avg{k2} = zeros(2,2);
    for k1 = 1:M
        psd_avg{k2} = psd_avg{k2} + psd_out{k2,k1};
    end
    psd_avg{k2} = psd_avg{k2}/M;
    %     disp(k2/length(f));
end


psd_plot = zeros(size(psd_avg,1),1);
for k3 = 1:size(psd_plot,1)
    psd_plot(k3) = real(psd_avg{k3}(1,1));
end
end