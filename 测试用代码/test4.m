%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2022-07-13 10:56:26
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2022-08-22 20:34:30
%FilePath: \NonlinearScanlan\test4.m
%Description: 模态分析
%
%Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% load kkmm.mat
% 
% % 特征值分析，即计算频率Freq和振型Phi，calmodes数字代表求解的阶数，eigs中参数SM表示从较小的特征值开始求解
% % Eigenvalue analysis, that is to calculate the frequency Freq and mode shape Phi, the calmodes number represents the order of the solution, and the parameter SM in eigs represents the solution from the smaller eigenvalue.
% calmodes = 5; %考虑模态数 Consider the number of modes
% [eig_vec, eig_val] = eigs(KK, MM, calmodes, 'SM');
% [nfdof, nfdof] = size(eig_vec);
% [omeg, w_order] = sort(sqrt(diag(eig_val)));
% mode_vec = eig_vec(:, w_order);
% Freq = omeg / (2 * pi);

load mode.mat
line = [1 2 3 4 5];
mode_vec(:,4)=mode_vec(:,4)*-1;
mode_vec(:,5)=mode_vec(:,5)*-1;
plot(line,mode_vec(1:5,:))
legend('1','2','3','4','5')


hold on

test=[-0.000188064841850800	0.00660936675881170	0.000207274293981572	0.00201556667806285	0.00351694075310958];
test=test/max(abs(test))*-1;
plot(line,test)
% KKK=KK(1:2,1:2);
% KKK=KKK*4;
% KKK(1,1)=KKK(1,1)/4;
% 
% MMM=MM(1:2,1:2);
% MMM(2,2)=MMM(2,2)*4;
% 
% calmodes = 2; %考虑模态数 Consider the number of modes
% [eig_vec, eig_val] = eigs(KKK, MMM, calmodes, 'SM');
% [nfdof, nfdof] = size(eig_vec);
% 
% [omeg, w_order] = sort(sqrt(diag(eig_val)));
% mode_vec = eig_vec(:, w_order);
% Freq = omeg / (2 * pi);