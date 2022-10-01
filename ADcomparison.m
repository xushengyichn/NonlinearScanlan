%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-10-01 10:05:10
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-10-01 11:10:51
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

figure
plot(Amplitude_exp*D_exp,zeta_exp)


function    [Zeta]=polynomial_zeta_exp(Amplitude,a1,a2,a3,a4,a5,rho,U,D,omega0,m,L)
    Zeta= -rho .* U .* D .* (a1 + 4 .* a2 .* Amplitude ./ 3 ./ pi + a3 .* Amplitude.^2/4 + 8 .* a4 .* Amplitude.^3/15 / pi + a5 .* Amplitude.^4/8)*L ./ 2 ./ omega0 ./ m;
end