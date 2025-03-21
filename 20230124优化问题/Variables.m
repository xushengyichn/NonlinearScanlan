%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 15:24:39
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-29 21:04:56
%FilePath: \20230124优化问题\Variables.m
%Description: 本脚本功能为创建所需优化的变量
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
addpath("../函数/")
% 输入预先设定的参数
% Preset parameters

modeinfo_all=importdata('modeinfo_all.mat');
mode_re=modeinfo_all.mode_re;
mode_re1=mode_re(:,1);
nodegap=modeinfo_all.nodegap;
nodedeck=0:1:660;
phideck = interp1(nodegap, mode_re1, nodedeck, 'linear', 'extrap'); %插值以后任意点的振型

[~,loc_seq]=find(abs(phideck)>0.5)
phideck_large=phideck(loc_seq);
nodedeck_large=nodedeck(loc_seq);

figure
plot(nodedeck_large,phideck_large)

% 
% number_of_modes_to_control = [1, 2, 3, 4, 5, 6]; % 需要控制的模态数 The number of modes to be controlled
% number_of_modes_to_consider = 10; % 考虑的总模态数 The total number of modes considered
% number_of_tmds = 3; % 考虑的总TMD数 The total number of TMDs
% total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
% mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
% total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
% modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
% t_length = 100; % 设定计算时间长度
% 
% % 输入初始参数
% % Input initial parameters
% TMDs_mass = [total_tmd_mass / number_of_tmds, total_tmd_mass / number_of_tmds]; % TMD质量 The mass of TMDs
% TMDs_mass = [TMDs_mass, total_tmd_mass - sum(TMDs_mass)]; % 最后一个TMD的质量 The mass of the last TMD
% TMDs_frequency = [0.8, 1, 1.2]; % TMD频率 The frequency of TMDs
% TMDs_damping_ratio = [0.08, 0.09, 0.07]; % TMD阻尼 The damping of TMDs
% TMDs_location = [55, 155, 255]; % TMD位置 The location of TMDs

% [result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
