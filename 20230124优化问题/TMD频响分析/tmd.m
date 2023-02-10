%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2023-02-10 00:04:28
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2023-02-10 00:04:40
%FilePath: \NonlinearScanlan\20230124优化问题\TMD频响分析\tmd.m
%Description: 绘制单个TMD的频响函数
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;



R = sqrt(((-g ^ 6 + (1 + (-4 * mu * xi__2 ^ 2 - 4 * xi__2 ^ 2 + mu + 2) * f ^ 2) * g ^ 4 - ((mu + 1) * f ^ 2 - 4 * xi__2 ^ 2 + 2) * f ^ 2 * g ^ 2 + f ^ 4) ^ 2) / ((g ^ 8) + (0.4e1 * (mu + 1) * ((mu * xi__2 ^ 2) + (xi__2 ^ 2) - 0.1e1 / 0.2e1) * (f ^ 2) + (8 * f * mu * xi__1 * xi__2) + (4 * xi__1 ^ 2) - 0.2e1) * (g ^ 6) + ((1 + (mu + 1) ^ 2 * f ^ 4 + (16 * xi__1 ^ 2 * xi__2 ^ 2 - 8 * mu * xi__2 ^ 2 - 8 * xi__1 ^ 2 - 8 * xi__2 ^ 2 + 2 * mu + 4) * f ^ 2) * g ^ 4) - (2 * (f ^ 2 * (-2 * xi__1 ^ 2 + mu + 1) - 2 * xi__2 ^ 2 + 1) * f ^ 2 * g ^ 2) + (f ^ 4)) ^ 2 + 0.4e1 * (((xi__2 * f * mu + xi__1) * g ^ 4) + 0.4e1 * ((xi__2 ^ 2) - 0.1e1 / 0.2e1) * xi__1 * (f ^ 2) * (g ^ 2) + (f ^ 4 * xi__1)) ^ 2 * (g ^ 2) / ((g ^ 8) + (0.4e1 * (mu + 1) * (-0.1e1 / 0.2e1 + ((mu + 1) * xi__2 ^ 2)) * (f ^ 2) + (8 * f * mu * xi__1 * xi__2) + (4 * xi__1 ^ 2) - 0.2e1) * (g ^ 6) + ((1 + (mu + 1) ^ 2 * f ^ 4 + ((16 * xi__1 ^ 2 - 8 * mu - 8) * xi__2 ^ 2 - 8 * xi__1 ^ 2 + 2 * mu + 4) * f ^ 2) * g ^ 4) - (2 * (f ^ 2 * (-2 * xi__1 ^ 2 + mu + 1) - 2 * xi__2 ^ 2 + 1) * f ^ 2 * g ^ 2) + (f ^ 4)) ^ 2);

