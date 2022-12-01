%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Shengyi xushengyichn@outlook.com
%Date: 2022-11-28 17:39:07
%LastEditors: Shengyi xushengyichn@outlook.com
%LastEditTime: 2022-11-28 17:40:21
%FilePath: \NonlinearScanlan\20221114一阶模态二TMD结果穷举\data_analysis.m
%Description: 分析数据
%
%Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;

% 读取数据
data = importdata('nmodes_onetmd_results_loc.mat');

loc=data(1:661,1);
dis=data(1:661,3);
figure
plot(loc,dis)

loc2=data(661*1+1:661*2,1);
dis2=data(661*1+1:661*2,3);
hold on
plot(loc2,dis2)

loc3=data(661*2+1:661*3,1);
dis3=data(661*2+1:661*3,3);
hold on
plot(loc3,dis3)

loc4=data(661*3+1:661*4,1);
dis4=data(661*3+1:661*4,3);
hold on
plot(loc4,dis4)

loc5=data(661*4+1:661*5,1);
dis5=data(661*4+1:661*5,3);
hold on
plot(loc5,dis5)
