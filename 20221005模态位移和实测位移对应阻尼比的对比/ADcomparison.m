%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-10-01 10:05:10
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-03 18:10:15
%FilePath: \NonlinearScanlan\ADcomparison.m
%Description: aerodynamic damping ratio comparison between real displacement and mode displacement
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear 
close all


% 导入试验数据
my_table = load('SZTD110_logfile.mat');
my_table = my_table.my_table;
Name = my_table.casename;
isexist = find(Name == 'SZTD-110-case2-22.3-fasan-2401');
a1 = my_table.up_parameter_a1(isexist); 
a2 = my_table.up_parameter_a2(isexist);
a3 = my_table.up_parameter_a3(isexist);
a4 = my_table.up_parameter_a4(isexist);
a5 = my_table.up_parameter_a5(isexist);
%% 试验阻尼比计算与实桥计算对比
m_exp= 80;
U_exp=6.21;
D_exp=0.667;
omega0_exp=32.8990;
Amplitude_exp=0:0.0001:0.01;
rho = 1.225;
L=3.6;
[zeta_exp]=polynomial_zeta_exp(Amplitude_exp,a1,a2,a3,a4,a5,rho,U_exp,D_exp,omega0_exp,m_exp,L);

% figure
% plot(Amplitude_exp,zeta_exp)
% 
% figure
% plot(Amplitude_exp*D_exp,zeta_exp)

%% mode shape of a pin-pin beam
dl=0.01
z=[0.01:dl:1];z=z';
omega0=omega0_exp;
EI=1;
rhoL=m_exp/max(z);
mode_number=1;
beta=mode_number*pi/max(z);
mode(:,1)=sin(beta.*z);
 

m_modal_temp=sum(rhoL*mode(:,1).^2)*dl;

mode_re=mode(:,1)/sqrt(m_modal_temp);
% mode_re=mode(:,1);

m_modal=sum(rhoL*mode_re(:,1).^2)*dl;
m_integral_1=sum(rhoL*mode_re(:,1))*dl;
m_integral_2=sum(rhoL*mode_re(:,1).^2)*dl;
m_integral_3=sum(rhoL*mode_re(:,1).^3)*dl;
m_integral_4=sum(rhoL*mode_re(:,1).^4)*dl;
m_integral_5=sum(rhoL*mode_re(:,1).^5)*dl;
m_integral_6=sum(rhoL*mode_re(:,1).^6)*dl;
integral_1=sum(abs(mode_re(:,1)))*dl;
integral_2=sum(mode_re(:,1).^2)*dl;
integral_3=sum(abs(mode_re(:,1)).^3)*dl;
integral_4=sum(mode_re(:,1).^4)*dl;
integral_5=sum(abs(mode_re(:,1)).^5)*dl;
integral_6=sum(mode_re(:,1).^6)*dl;
integral_7=sum(abs(mode_re(:,1)).^7)*dl;

% figure
% plot(z,mode_re)
U=U_exp;
D=D_exp;
Amplitude1=Amplitude_exp/D;%试验无量纲振幅
Amplitude2=Amplitude_exp/max(mode_re)/D;%模态坐标无量纲振幅
m=m_modal;

%模态坐标下气动力系数换算
d1=m_modal*a1/m_exp;
d2=m_modal*a2/m_exp*max(mode_re);
d3=m_modal*a3/m_exp*max(mode_re)^2;
d4=m_modal*a4/m_exp*max(mode_re)^3;
d5=m_modal*a5/m_exp*max(mode_re)^4;



% zeta_exp= -rho .* U_exp .* D_exp .* (a1 +4 .* a2 .* Amplitude*max(mode_re) ./ 3 ./ pi+a3.*(Amplitude*max(mode_re)).^2/4+8*a4.*(Amplitude*max(mode_re)).^3/15/pi+a5.*(Amplitude*max(mode_re)).^4/8) ./ 2 ./ omega0_exp ./ m_exp;

% Zeta=-rho.*U.*D.*(d1+4*d2.*Amplitude/3/pi+d3.*Amplitude.^2/4+8*d4.*Amplitude.^3/15/pi+d5.*Amplitude.^4/8)./ 2 ./ omega0 ./m_modal;


zeta_exp= -rho .* U_exp .* D_exp .* (a1 +4 .* a2 .* Amplitude1 ./ 3 ./ pi+a3.*Amplitude1.^2/4+8.*a4.*(Amplitude1).^3/15/pi+a5.*Amplitude1.^4/8) ./ 2 ./ omega0_exp ./ m_exp

Zeta=-rho.*U.*D.*(d1+4*d2.*Amplitude2/3/pi+d3.*Amplitude2.^2/4+8.*d4.*(Amplitude2).^3/15/pi+d5.*Amplitude2.^4/8)./ 2 ./ omega0 ./m_modal


figure

plot(Amplitude_exp,Zeta)
xlim([0 0.006]);
hold on
plot(Amplitude_exp,zeta_exp)

legend("modal coordinates","rigid model")

% %%%mode shape of a 6-span beam
% lamda=[3.141592654	3.260534994	3.55640846	3.926602312	4.297529693	4.601417878];
% z_virtue=[0.01:0.01:1];z_virtue=z_virtue';
% z=[0:0.01:6];z=z';

% for ii=1:6
%     A(1)=1;B(1)=0;
%     mode_virtue(:,1)=A(1)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(1)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
    
%     for jj=2:6
%         A(jj)=2*A(jj-1)*(sinh(lamda(ii))*cos(lamda(ii))-sin(lamda(ii))*cosh(lamda(ii)))/(sinh(lamda(ii))-sin(lamda(ii)))-B(jj-1);
%         B(jj)=A(jj-1);
%         mode_virtue(:,jj)=A(jj)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(jj)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
%     end
%     mode(:,ii)=[mode_virtue(:,1);mode_virtue(:,2);mode_virtue(:,3);mode_virtue(:,4);mode_virtue(:,5);mode_virtue(:,6);];
% %     mode(:,ii)=[mode_virtue(:,1)];
%     mode(:,ii)=mode(:,ii)/max(mode(:,ii));
% end
% figure
% plot(mode(:,1))

function    [Zeta]=polynomial_zeta_exp(Amplitude,a1,a2,a3,a4,a5,rho,U,D,omega0,m,L)
    Zeta= -rho .* U .* D .* (a1 + 4 .* a2 .* Amplitude ./ 3 ./ pi + a3 .* Amplitude.^2/4 + 8 .* a4 .* Amplitude.^3/15 / pi + a5 .* Amplitude.^4/8)*L ./ 2 ./ omega0 ./ m;
end