%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-10-01 10:05:10
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-10-01 15:32:42
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
% 
% figure
% plot(Amplitude_exp*D_exp,zeta_exp)


%%%mode shape of a 6-span beam
lamda=[3.141592654	3.260534994	3.55640846	3.926602312	4.297529693	4.601417878];
z_virtue=[0.01:0.01:1];z_virtue=z_virtue';
z=[0:0.01:6];z=z';

for ii=1:6
    A(1)=1;B(1)=0;
    mode_virtue(:,1)=A(1)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(1)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
    
    for jj=2:6
        A(jj)=2*A(jj-1)*(sinh(lamda(ii))*cos(lamda(ii))-sin(lamda(ii))*cosh(lamda(ii)))/(sinh(lamda(ii))-sin(lamda(ii)))-B(jj-1);
        B(jj)=A(jj-1);
        mode_virtue(:,jj)=A(jj)*(sinh(lamda(ii))*sin(lamda(ii)*z_virtue)-sin(lamda(ii))*sinh(lamda(ii)*z_virtue))+B(jj)*(sinh(lamda(ii))*sin(lamda(ii)*(1-z_virtue))-sin(lamda(ii))*sinh(lamda(ii)*(1-z_virtue)));
    end
    mode(:,ii)=[mode_virtue(:,1);mode_virtue(:,2);mode_virtue(:,3);mode_virtue(:,4);mode_virtue(:,5);mode_virtue(:,6);];
%     mode(:,ii)=[mode_virtue(:,1)];
    mode(:,ii)=mode(:,ii)/max(mode(:,ii));
end
figure
plot(mode(:,1))

function    [Zeta]=polynomial_zeta_exp(Amplitude,a1,a2,a3,a4,a5,rho,U,D,omega0,m,L)
    Zeta= -rho .* U .* D .* (a1 + 4 .* a2 .* Amplitude ./ 3 ./ pi + a3 .* Amplitude.^2/4 + 8 .* a4 .* Amplitude.^3/15 / pi + a5 .* Amplitude.^4/8)*L ./ 2 ./ omega0 ./ m;
end