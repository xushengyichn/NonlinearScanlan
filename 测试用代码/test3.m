clc
clear
close all

frequencyOffset=(-2:0.1:2)/100;

% for k1 = 1:length(massoftmd)
displacement=zeros(length(frequencyOffset),1);
Frequency=zeros(length(frequencyOffset),1);
parfor k1 = 1:length(frequencyOffset)
    [out]=test_function(frequencyOffset(k1));
    dis=out(:,2);
    disend=dis(end/2:end);
    displacement(k1,:)=[max(disend)];
    [~,a,b]=fft_transfer(256,dis)
    [~,index]=max(b)
    Frequency(k1,:)=a(index);
end

figure
plot(frequencyOffset,displacement(:,1))
figure
plot(frequencyOffset,Frequency(:,1))
% hold on
% plot(massoftmd,displacement(:,2))