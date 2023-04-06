%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-16 11:37:30
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-04-06 10:19:25
%FilePath: \NonlinearScanlan\20221007双幅桥多项式模型参数识别(增加振幅上下限限制点，避免分段函数)\twindeck_VIV.m
%Description: 提取双幅桥振动响应，并使用多项式模型进行参数识别

%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
% addpath('..\函数\')
%% 各种参数设置
recorddata = 1; % 是否需要记录数据
fitdecide =0; % 是否需要拟合振幅包络线
fname = ['SZTD-110-case2-22.3-fasan-2401']; %分析工况对应文件名
fname_up_dltx= ['SZTD-110-case2-22.3-dltx2'];%上游动力特性文件名
fname_down_dltx= ['SZTD-110-case2-22.3-dltx4'];%下游动力特性文件名
L = 3.6; % 节段模型长度
m_up=80/L;%上游梁单位长度质量
m_down=80/L;%下游梁单位长度质量
windspeeds=importdata('windspeeds.mat');
windspeed=windspeeds.Windspeed(find(windspeeds.casename==fname));
disp("该工况的风速为"+windspeed+"m/s")

my_table=importdata("SZTD110_logfile.mat");

casenumber=find(my_table.casename==fname);


% 结构参数（上游）
% Structural parameters （windward side）
D = 0.667; % deck depth
m = 80/L; % mass of the segment model per unit length
F0 = my_table.up_Fre_vibration; % Frequency without wind
omega0 = 2 * pi * F0; % Circular frequency without wind
rho = 1.225; % density of the air
% U = 6.21; % wind speed
Name = my_table.casename;
isexist = find(Name == fname);
if ~isempty(my_table.Windspeed(isexist))
    alertstr=(fname+"工况风速为："+num2str(my_table.Windspeed(isexist)));
    U=my_table.Windspeed(isexist);
else
        alertstr=("工况名可能输入错误");
end
% Zeta0 = 0.29/100; % damping ratio without wind



Zeta0 = 0.003;
up_a1 = my_table.up_parameter_a1(casenumber);
up_a2 = my_table.up_parameter_a2(casenumber);
up_a3 = my_table.up_parameter_a3(casenumber);
up_a4 = my_table.up_parameter_a4(casenumber);
up_a5 = my_table.up_parameter_a5(casenumber);

up_a = [up_a1 up_a2 up_a3 up_a4 up_a5];
up_H4= my_table.up_parameter_H4(casenumber);


Fre = my_table.up_Fre_vibration(casenumber);
Mass = m;


up_told = (0:0.01:100)';

fsup_a = [up_a1 up_a2 up_a3 up_a4 up_a5];
up_P = zeros(1, size(up_told, 1));


up_u0 = 0.001;
up_udot0 = 0;

up_out = polynomial_NB_adstiff(Fre, Mass, Zeta0, rho, D, U, up_a, up_H4, up_told, up_P, up_u0, up_udot0);



% 结构参数（下游）
% Structural parameters （leeward side）
F0 = my_table.down_Fre_vibration; % Frequency without wind
omega0 = 2 * pi * F0; % Circular frequency without wind



Zeta0 = 0.003;
down_a1 = my_table.down_parameter_a1(casenumber);
down_a2 = my_table.down_parameter_a2(casenumber);
down_a3 = my_table.down_parameter_a3(casenumber);
down_a4 = my_table.down_parameter_a4(casenumber);
down_a5 = my_table.down_parameter_a5(casenumber);

down_a = [down_a1 down_a2 down_a3 down_a4 down_a5];
down_H4=my_table.down_parameter_H4(casenumber);

Fre = my_table.down_Fre_vibration(casenumber);
Mass = m;
down_u0 = 0.001;
down_udot0 = 0;



down_told = (0:0.01:100)';
fsdown_a = [down_a1 down_a2 down_a3 down_a4 down_a5];
down_P = zeros(1, size(down_told, 1));



down_out = polynomial_NB_adstiff(Fre, Mass, Zeta0, rho, D, U, down_a,down_H4, down_told, down_P, down_u0, down_udot0);



figure
plot(up_out(:, 1), up_out(:, 2))
hold on
plot(down_out(:, 1), down_out(:, 2))
legend("windward girder","leeward girder")
title("振动时程重构")
