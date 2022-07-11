clc
clear
close all

massoftmd=10:10:250;
massoftmd=massoftmd';
massoftmd=massoftmd/1000;

% for k1 = 1:length(massoftmd)
displacement=zeros(length(massoftmd),2);
parfor k1 = 1:length(massoftmd)
    [out,out1]=test_function(massoftmd(k1));
    dis=out(:,2);
    disend=dis(end/2:end);
    
    dis1=out1(:,2);
    disend1=dis1(end/2:end);
    displacement(k1,:)=[max(disend) max(disend1)];
end

figure
plot(massoftmd,displacement(:,1))
% hold on
% plot(massoftmd,displacement(:,2))