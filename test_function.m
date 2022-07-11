%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-27 16:21:36
%LastEditors: Shengyi Xu xushengyichn@outlook.com
%LastEditTime: 2022-07-11 14:33:50
%FilePath: \NonlinearScanlan\Polynomial_withTMD_singledegree.m
%Description: 本代码是用于求解多项式模型下，单自由度模型安装TMD后的响应计算。
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clc
% clear
% This is a function
function [out,out1] = test_function(massoftmd)

% Nonlinear Newmark's Direct Integration Method with polynomial model
% (n = number of time steps)
% (ndof = number degrees of freedom)

% INPUT
% Fre      = Frequency of the system         =>[1]
% Mass     = Mass of the system              =>[1]
% Zeta0    = Damping ratio of the system     =>[1]
% rho      = Air density                     =>[1]
% D        = Reference length                =>[1]
% U        = Wind speed at a certain reduced frequency       =>[1]
% Y1k      = Y1 at a certain reduced frequency       =>[1]
% epsilonk = Epsilon at a certain reduced frequency       =>[1]
% Y2k      = Y2 at a certain reduced frequency       =>[1]
% Clk      = Cl at a certain reduced frequency       =>[1]
% t        = Time vector         =>[1,n]
% P        = load vs. time       =>[1,n]
% u0       = Initial displacements =>[1]
% udot0    = Initial velocity =>[1]
% gam      = gamma (constant)
% beta     = beta  (constant)
% mtmd     = mass of tmd
% ctmd     = damping of tmd
% ktmd     = stiffness of tmd

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 => explicit central difference method
% beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
gamma = 1/2; % Factor in the Newmark algorithm
beta = 1/4; % Factor in the Newmark algorithm
%--------------------------------------------------------------------------
% 设置矩阵大小
matrixsize = 2; %calculate one degree of freedom
%--------------------------------------------------------------------------
% 导入试验数据
load('TMD_logfile.mat');
my_table_tmd=my_table;
load('SZTD110_logfile.mat');

%% 提取工况设置
fname = ['SZTD-110-case2-22.3-fasan-2401']; %记录文件名

chanum = 8; %记录通道数
filename = strsplit(fname, '-');
casenumber = cell2mat(filename(3));
spacing = str2double(cell2mat(filename(4)));
type = cell2mat(filename(5));
Rspeed = cell2mat(filename(6));

casenumber = str2double(casenumber(5:end));
Rspeed = Rspeed(1:end - 1);
Rspeed = str2double(Rspeed);

if strcmp(type, 'fasan')
    type = 1;
else

    if strcmp(type, 'shuaijian')
        type = 2;
    else
        error("实验数据文件名可能有误，请确认信号状态为衰减或发散，程序终止")
    end

end

% 判断工况风攻角
if casenumber <= 17 || and(casenumber >= 32, casenumber <= 35) || and(casenumber >= 40, casenumber <= 16) || casenumber == 55
    AOA = 3;
else

    if and(casenumber >= 18, casenumber <= 25) || and(casenumber >= 36, casenumber <= 37) || and(casenumber >= 47, casenumber <= 50) || casenumber == 54
        AOA = 0;
    else

        if and(casenumber >= 26, casenumber <= 31) || and(casenumber >= 38, casenumber <= 39) || and(casenumber >= 51, casenumber <= 53) || and(casenumber >= 56, casenumber <= 57)
            AOA = -3;
        end

    end

end

%导入试验数据
Name = my_table.casename;
isexist = find(Name == fname);
% up_index_start = my_table.up_start(isexist);
% up_index_end = my_table.up_end(isexist);
% down_index_start = my_table.down_start(isexist);
% down_index_end = my_table.down_end(isexist);
up_a1=my_table.up_parameter_a1(isexist);
up_a2=my_table.up_parameter_a2(isexist);
up_a3=my_table.up_parameter_a3(isexist);
up_a4=my_table.up_parameter_a4(isexist);
up_a5=my_table.up_parameter_a5(isexist);


up_a = [up_a1 up_a2 up_a3 up_a4 up_a5];

up_H4=my_table.up_parameter_H4(isexist);% 气动刚度
up_upperlimit= my_table.up_upperlimit(isexist); %除以特征长度D的无量纲振幅
up_lowerlimit= my_table.up_lowerlimit(isexist); %除以特征长度D的无量纲振幅
up_Fren_vibration_withwind=my_table.up_Fren_vibration_withwind(isexist); 


% 结构参数
% Structural parameters
D = 0.667; % deck depth
m = 80; % mass of the segment model
Mass = m;
% F0 = 5.276; % Frequency without wind
F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
Fre= F0;
% omega0 = 2 * pi * F0; % Circular frequency without wind
rho = 1.225; % density of the air
U = my_table.Windspeed(isexist); % density of the air
% U = 6.21; % wind speed
% Zeta0 = 0.3/100; % damping ratio without wind
Zeta0 =  my_table.up_dltx_zeta0(isexist); % damping ratio without wind

h = 1/256; % 时间步长 % Step size of the algorithm
up_t=0:h:60; % Time vector
up_told = up_t;
up_tt = up_told;
up_P = zeros(1, size(up_told, 1));
P = zeros(2, length(up_tt));


% up_u0 = [-6.3271e-04; 0];
% up_udot0 = [-0.0161; 0];

up_u0 = [-6.3271e-04; 0;0;0;0];
up_udot0 = [-0.0161; 0;0;0;0];
P = zeros(5, length(up_tt));
%% TMD 参数
% sel=[11 13 15 17];
sel=[11 11 11 11];
zetatmd = my_table_tmd.zeta(sel);
fretmd = my_table_tmd.fre(sel);
mtmd = ones(length(sel),1)*massoftmd;

%% Calculate the response

% nModes = 1;
% matrixsize=2;
% out = polynomial_NB_withTMD(Fre, Mass, Zeta0, rho, D, U, up_a, up_t,h, P,  up_u0, up_udot0,nModes,matrixsize);
up_u0 = [-6.3271e-04; 0];
up_udot0 = [-0.0161; 0];
P = zeros(2, length(up_tt));
nModes = 1;
matrixsize=2;
out1 = test4(Fre, Mass, Zeta0, rho, D, U, up_a, up_t,h, P,  up_u0, up_udot0,nModes,matrixsize,massoftmd*4);

up_u0 = [-6.3271e-04; 0;0;0;0];
up_udot0 = [-0.0161; 0;0;0;0];
P = zeros(5, length(up_tt));
nModes = 1;
matrixsize=5;
% out = polynomial_NB_withTMDs(Fre, Mass, Zeta0, rho, D, U, up_a, up_t,h, P, up_u0, up_udot0,nModes,mtmd,fretmd,zetatmd);

out = polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);

% nModes = 1;
% matrixsize=5;
% % out = polynomial_NB_withTMDs(Fre, Mass, Zeta0, rho, D, U, up_a, up_t,h, P, up_u0, up_udot0,nModes,mtmd,fretmd,zetatmd);
% 
% out = polynomial_NB_withTMDs_addstiff_withlimit(Fre, Mass, Zeta0, rho, D, U, up_a,up_H4, up_t,h, P, up_u0, up_udot0,up_upperlimit,up_lowerlimit,up_Fren_vibration_withwind,nModes,mtmd,fretmd,zetatmd);
% 
% subplot(1, 2, 1)
% plot(out(:, 1), out(:, 2))
% ylim([-0.01 0.01])
% title("main structure")
% xlabel('time(s)');
% ylabel('displacement(m)');
% subplot(1, 2, 2)
% plot(out(:, 1), out(:, 3))
% title("TMD")
% ylim([-0.01 0.01])

figure 
plot(out(:, 1), out(:, 2))
ylim([-0.01 0.01])
title("main structure")
hold on
load test.mat
plot(up_t,UP)
legend("calculated","windtunnel test")
fs=1/(up_t(2)-up_t(1));
[psd_avg, f1, psd_plot1] = fft_transfer(fs,UP);
[psd_avg, f2, psd_plot2] = fft_transfer(256,out(:, 2));
% figure
% plot(f1, psd_plot1)
% hold on 
% plot(f2, psd_plot2)
% legend("windtunnel test","calculated")
max(out(:, 2))
[a,b]=max(psd_plot1);
f1(b)

[a,b]=max(psd_plot2);
f2(b)
close all
end
