%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi Xu xushengyichn@outlook.com
%Date: 2022-06-27 16:21:36
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-09-14 10:19:09
%FilePath: \NonlinearScanlan\test_Polynomial_withTMD_singledegree.m
%Description: 本代码是用于求解多项式模型下，测试节段模型安装tmd后产生多阶模态，最后计算得到响应的频率问题。
%
%Copyright (c) 2022 by Shengyi Xu xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
clear


%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 => explicit central difference method
% beta = 1/4,   gamma = 1/2 => undamped trapezoidal rule (implicit)

%--------------------------------------------------------------------------
gamma = 1/2; % Factor in the Newmark algorithm
beta = 1/4; % Factor in the Newmark algorithm
%--------------------------------------------------------------------------
% 设置矩阵大小
matrixsize = 5; %calculate one degree of freedom
%--------------------------------------------------------------------------

% 结构参数
% Structural parameters
D = 0.667; % deck depth
L = 3.6; %length of the sectional model
m = 80; % mass of the segment model per unit length
Mass = m;
F0 = 0.83; % Frequency without wind
Fre= F0;
rho = 1.225; % density of the air

Zeta0 =  0.003; % damping ratio without wind

h = 1/256; % 时间步长 % Step size of the algorithm
t=0:h:60; % Time vector

p=sin(2*pi*F0*t);


%% TMD 参数
sel=[12 14 15 18];
% sel=[4 7 8 10];
% sel=[11 11 11 11];
% sel=[11];
zetatmd = [0.1 0.1 0.1 0.1];
fretmd = [0.83 0.84 0.85 0.86];
mtmd = ones(length(sel),1)*0.25;
disp("TMD质量："+num2str(mtmd));
disp("TMD阻尼系数："+num2str(zetatmd));
disp("TMD频率："+num2str(fretmd));
%% Calculate the response


u0 = [-0.086; -1;-1;-1;-1]*10e-3;
udot0 = [0; 0;0;0;0];

P = zeros(5, length(t));
P(1,:)=p;
nModes = 1;
matrixsize=5;



[MM,CC,KK]=CreateMatrixwithTMD(nModes,Mass,Zeta0,Fre,mtmd,zetatmd,fretmd);


[u udot u2dot]=NewmarkInt(t,MM,CC,KK,P,gamma,beta,u0,udot0);


figure
subplot(1, 2, 1)
plot(t, u(1,:))
ylim([-0.05 0.05])
title("main structure")
xlabel('time(s)');
ylabel('displacement(m)');
subplot(1, 2, 2)
plot(t, u(2,:))
title("TMD1")
ylim([-0.05 0.05])
