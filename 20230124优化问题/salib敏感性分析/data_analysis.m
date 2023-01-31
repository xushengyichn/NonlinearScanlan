%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-30 10:19:36
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-30 10:20:04
%FilePath: \20230124优化问题\salib敏感性分析\data_analysis.m
%Description: 将计算结果绘图
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% 读取数据
input = importdata('param_values_mode1.txt');
output = importdata('salib_mode1.txt');

% 绘制散点图
figure
scatter(input(:,2),output(:))