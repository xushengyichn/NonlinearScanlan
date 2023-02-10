clc
clear
close all

m1 = 100;
xi1 = 0.003;
f1=1;
w1=2*pi*f1;
c1=2*m1*w1*xi1;
k1=m1*w1^2;

m2 = 1;
xi2= 0.003;
f2 = 1;
w2 = 2*pi*f2;
c2=2*m2*w2*xi2;
k2=m2*w2^2;

f_F=0.001:0.01:2;


F=1;
for i=1:length(f_F)
    w_F(i)=2*pi*f_F(i);
    w=w_F(i)
%     X(i)=F/(m1*w_F(i)^2+(c1+c2)*w_F(i)+(k1+k2)-(k2+c2*w_F(i))^2/(m2*w_F(i)^2+c2*w_F(i)+k2));
%     Xs=F/k1;
%     result(i)=X(i)/Xs;
    G(i)=(-m2*w^2-c2*w+k2)*k1/(-m1*w^2-(c1+c2)*w+k1+k2)*(-m2*w^2-c2*w+k2)-(k2-c2*w)^2;
end
plot(f_F,G)
